#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
AWS-related methods.
"""

import os.path as op
import sys
import fnmatch

import boto3
import json
import logging
import time
from multiprocessing import Pool

from jcvi.formats.base import BaseFile, SetFile, timestamp
from jcvi.apps.base import OptionParser, ActionDispatcher, datafile, popen, sh


class InstanceSkeleton(BaseFile):

    def __init__(self, filename=datafile("instance.json")):
        super(InstanceSkeleton, self).__init__(filename)
        self.spec = json.load(open(filename))

    @property
    def launch_spec(self):
        return self.spec["LaunchSpec"]

    @property
    def instance_id(self):
        return self.spec["InstanceId"]

    @property
    def volumes(self):
        return self.spec["Volumes"]

    @property
    def image_id(self):
        return self.spec["LaunchSpec"]["ImageId"]

    def save(self):
        fw = open(self.filename, "w")
        json.dump(self.spec, fw, indent=4)
        fw.close()

    def save_instance_id(self, instance_id):
        self.spec["InstanceId"] = instance_id
        self.save()

    def save_image_id(self, image_id):
        self.spec["LaunchSpec"]["ImageId"] = image_id
        self.save()


def main():

    actions = (
        ('cp', 'copy files with support for wildcards'),
        ('ls', 'list files with support for wildcards'),
        ('rm', 'remove files with support for wildcards'),
        ('role', 'change aws role'),
        ('launch', 'launch ec2 instance'),
        ('stop', 'stop ec2 instance'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def launch(args):
    """
    %prog launch

    Launch ec2 instance through command line.
    """
    p = OptionParser(launch.__doc__)
    p.add_option("--profile", default="mvrad-datasci-role", help="Profile name")
    p.add_option("--price", default=3.0, type=float, help="Spot price")
    opts, args = p.parse_args(args)

    if len(args) != 0:
        sys.exit(not p.print_help())

    role(["205134639408", "htang", "114692162163", "mvrad-datasci-role"])
    session = boto3.Session(profile_name=opts.profile)
    client = session.client('ec2')
    s = InstanceSkeleton()

    # Make sure the instance id is empty
    instance_id = s.instance_id
    if len(instance_id) > 0:
        logging.error("Instance exists {}".format(instance_id))
        sys.exit(1)

    # Launch spot instance
    launch_spec = s.launch_spec
    response = client.request_spot_instances(
        SpotPrice=str(opts.price),
        InstanceCount=1,
        Type="one-time",
        AvailabilityZoneGroup="us-west-2b",
        LaunchSpecification=launch_spec
    )

    request_id = response["SpotInstanceRequests"][0]["SpotInstanceRequestId"]
    print >> sys.stderr, "Request id {}".format(request_id)

    instance_id = ""
    while not instance_id:
        response = client.describe_spot_instance_requests(
            SpotInstanceRequestIds=[request_id]
        )
        if "InstanceId" in response["SpotInstanceRequests"][0]:
            instance_id = response["SpotInstanceRequests"][0]["InstanceId"]
        else:
            logging.debug("Waiting to be fulfilled ...")
            time.sleep(10)

    print >> sys.stderr, "Instance id {}".format(instance_id)

    status = ""
    while status != "running":
        logging.debug("Waiting instance to run ...")
        time.sleep(3)
        response = client.describe_instance_status(InstanceIds=[instance_id])
        if len(response["InstanceStatuses"]) > 0:
            status = response["InstanceStatuses"][0]["InstanceState"]["Name"]

    # Tagging
    response = client.create_tags(
        Resources=[instance_id],
        Tags=[{"Key": k, "Value": v} for k, v in { \
                    "Name": "htang-lx-spot",
                    "owner": "htang",
                    "project": "mv-bioinformatics"
                }.items()]
    )

    # Attach working volumes
    volumes = s.volumes
    for volume in volumes:
        response = client.attach_volume(
            VolumeId=volume["VolumeId"],
            InstanceId=instance_id,
            Device=volume["Device"]
        )

    # Save instance id and ip
    s.save_instance_id(instance_id)
    response = client.describe_instances(InstanceIds=[instance_id])
    ip_address = response["Reservations"][0]["Instances"][0]["PrivateIpAddress"]
    print >> sys.stderr, "IP address {}".format(ip_address)


def stop(args):
    """
    %prog stop

    Stop EC2 instance.
    """
    p = OptionParser(stop.__doc__)
    p.add_option("--profile", default="mvrad-datasci-role", help="Profile name")
    opts, args = p.parse_args(args)

    if len(args) != 0:
        sys.exit(not p.print_help())

    role(["205134639408", "htang", "114692162163", "mvrad-datasci-role"])
    session = boto3.Session(profile_name=opts.profile)
    client = session.client('ec2')
    s = InstanceSkeleton()

    # Create image
    instance_id = s.instance_id
    block_device_mappings = []
    for volume in s.volumes:
        block_device_mappings.append(
            {
                "DeviceName": volume["Device"],
                "NoDevice": ""
            }
        )

    new_image_name = "htang-dev-{}-{}".format(timestamp(), int(time.time()))
    response = client.create_image(
        InstanceId=instance_id,
        Name=new_image_name,
        BlockDeviceMappings=block_device_mappings
    )
    print >> sys.stderr, response
    new_image_id = response["ImageId"]

    image_status = ""
    while image_status != "available":
        logging.debug("Waiting for image to be ready")
        time.sleep(10)
        response = client.describe_images(ImageIds=[new_image_id])
        image_status = response["Images"][0]["State"]

    # Delete old image, snapshot and shut down instance
    old_image_id = s.image_id
    response = client.describe_images(ImageIds=[old_image_id])
    old_snapshot_id = response["Images"][0]["BlockDeviceMappings"][0]["Ebs"]["SnapshotId"]
    response = client.deregister_image(ImageId=old_image_id)
    print >> sys.stderr, response
    response = client.delete_snapshot(SnapshotId=old_snapshot_id)
    print >> sys.stderr, response
    response = client.terminate_instances(InstanceIds=[instance_id])
    print >> sys.stderr, response

    # Save new image id
    s.save_image_id(new_image_id)
    s.save_instance_id("")


def glob_s3(store, keys=None, recursive=False):
    store, cards = store.rsplit("/", 1)
    contents = ls_s3(store, recursive=recursive)
    if keys:
        filtered = [x for x in contents if op.basename(x).split(".")[0] in keys]
    else:
        filtered = fnmatch.filter(contents, cards)

    if recursive:
        store = "s3://" + store.replace("s3://", "").split("/")[0]

    filtered = ["/".join((store, x)) for x in filtered]

    return filtered


def rm_s3(store):
    cmd = "aws s3 rm {}".format(store)
    sh(cmd)


def rm(args):
    """
    %prog rm "s3://hli-mv-data-science/htang/str/*.csv"

    Remove a bunch of files.
    """
    p = OptionParser(rm.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    store, = args
    contents = glob_s3(store)
    for c in contents:
        rm_s3(c)


def worker(work):
    c, target, force = work
    if force or not op.exists(target):
        pull_from_s3(c, target)


def cp(args):
    """
    %prog cp "s3://hli-mv-data-science/htang/str/*.csv" .

    Copy files to folder. Accepts list of s3 addresses as input.
    """
    p = OptionParser(cp.__doc__)
    p.add_option("--force", default=False,
                 action="store_true", help="Force overwrite if exists")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    store, folder = args
    force = opts.force
    cpus = opts.cpus
    if op.exists(store):
        contents = [x.strip().split(",") for x in open(store)]
    else:
        contents = glob_s3(store)

    tasks = []
    for c in contents:
        if isinstance(c, basestring):
            oc = op.basename(c)
            tc = op.join(folder, oc)
        else:
            c, tc = c
        tasks.append((c, tc, force))

    worker_pool = Pool(cpus)
    worker_pool.map(worker, tasks)
    worker_pool.close()
    worker_pool.join()


def ls(args):
    """
    %prog ls "s3://hli-mv-data-science/htang/str/*.vcf.gz"

    List files with support for wildcards.
    """
    p = OptionParser(ls.__doc__)
    p.add_option("--keys", help="List of keys to include")
    p.add_option("--recursive", default=False, action="store_true",
                 help="Recursive search")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    store, = args
    keys = opts.keys
    if keys:
        keys = SetFile(keys)
    print "\n".join(glob_s3(store, keys=keys, recursive=opts.recursive))


def s3ify(address):
    if not address.startswith("s3://"):
        address = "s3://" + address.lstrip("/")
    return address


def push_to_s3(s3_store, obj_name):
    cmd = "sync" if op.isdir(obj_name) else "cp"
    s3address = "{0}/{1}".format(s3_store, obj_name)
    s3address = s3ify(s3address)
    cmd = "aws s3 {0} {1} {2} --sse".format(cmd, obj_name, s3address)
    sh(cmd)
    return s3address


def pull_from_s3(s3_store, file_name=None, overwrite=True):
    is_dir = s3_store.endswith("/")
    if is_dir:
        s3_store = s3_store.rstrip("/")
    file_name = file_name or s3_store.split("/")[-1]
    if not op.exists(file_name):
        s3_store = s3ify(s3_store)
        if overwrite or (not op.exists(file_name)):
            cmd = "aws s3 cp {0} {1} --sse".format(s3_store, file_name)
            if is_dir:
                cmd += " --recursive"
            sh(cmd)
    return op.abspath(file_name)


def sync_from_s3(s3_store, target_dir=None):
    s3_store = s3_store.rstrip("/")
    s3_store = s3ify(s3_store)
    if target_dir is None:
        target_dir = op.basename(s3_store)
    cmd = "aws s3 sync {}/ {}/".format(s3_store, target_dir)
    sh(cmd)
    return target_dir


def ls_s3(s3_store_obj_name, recursive=False):
    s3_store_obj_name = s3ify(s3_store_obj_name)
    cmd = "aws s3 ls {0}/".format(s3_store_obj_name)
    if recursive:
        cmd += " --recursive"
    contents = []
    for row in popen(cmd):
        contents.append(row.split()[-1])
    return contents


def check_exists_s3(s3_store_obj_name):
    s3_store_obj_name = s3ify(s3_store_obj_name)
    cmd = "aws s3 ls {0} | wc -l".format(s3_store_obj_name)
    counts = int(popen(cmd).read())
    return counts != 0


def aws_configure(profile, key, value):
    sh('aws configure set profile.{0}.{1} {2}'.format(profile, key, value))


def role(args):
    """
    %prog role src_acct src_username dst_acct dst_role

    Change aws role. For example:
    %prog role 205134639408 htang 114692162163 mvrad-datasci-role
    """
    p = OptionParser(role.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 1 and args[0] == "htang":
        args = "205134639408 htang 114692162163 mvrad-datasci-role".split()

    if len(args) != 4:
        sys.exit(not p.print_help())

    src_acct, src_username, dst_acct, dst_role = args
    region = 'us-west-2'
    mfa_token = raw_input('Enter MFA Token [ENTER]: ')

    boto3.setup_default_session(profile_name='default')

    client = boto3.client('sts')
    response = client.assume_role(RoleArn="arn:aws:iam::" + dst_acct + ":role/" + dst_role,
                                  RoleSessionName=dst_role,
                                  SerialNumber="arn:aws:iam::" + src_acct + ":mfa/" + src_username,
                                  TokenCode=mfa_token)

    creds = response['Credentials']

    aws_configure(dst_role, 'aws_access_key_id', creds['AccessKeyId'])
    aws_configure(dst_role, 'aws_secret_access_key', creds['SecretAccessKey'])
    aws_configure(dst_role, 'aws_session_token', creds['SessionToken'])
    aws_configure(dst_role, 'region', region)

    print dst_role


if __name__ == '__main__':
    main()
