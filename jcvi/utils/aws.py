#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
AWS-related methods.
"""
import os
import os.path as op
import sys
import fnmatch
import boto3
import json
import logging
import time
import getpass

from configparser import NoOptionError, NoSectionError
from datetime import datetime
from multiprocessing import Pool
from botocore.exceptions import ClientError, ParamValidationError

from jcvi.formats.base import BaseFile, SetFile, timestamp
from jcvi.utils.console import console
from jcvi.apps.base import (
    OptionParser,
    ActionDispatcher,
    datafile,
    get_config,
    popen,
    sh,
)

AWS_CREDS_PATH = "%s/.aws/credentials" % (op.expanduser("~"),)


class InstanceSkeleton(BaseFile):
    def __init__(self, filename=datafile("instance.json")):
        super(InstanceSkeleton, self).__init__(filename)
        self.spec = json.load(open(filename))

    @property
    def launch_spec(self):
        return self.spec["LaunchSpec"]

    @property
    def instance_id(self):
        return self.spec["InstanceId"].strip()

    @property
    def private_ip_address(self):
        return self.spec["PrivateIpAddress"]

    @property
    def availability_zone(self):
        return self.spec["AvailabilityZone"]

    @property
    def volumes(self):
        return self.spec["Volumes"]

    @property
    def block_device_mappings(self):
        return self.launch_spec["BlockDeviceMappings"]

    @property
    def ebs_optimized(self):
        return self.launch_spec["EbsOptimized"]

    @property
    def image_id(self):
        return self.launch_spec["ImageId"]

    @property
    def instance_type(self):
        return self.launch_spec["InstanceType"]

    @property
    def key_name(self):
        return self.launch_spec["KeyName"]

    @property
    def security_group_ids(self):
        return self.launch_spec["SecurityGroupIds"]

    @property
    def subnet_id(self):
        return self.launch_spec["SubnetId"]

    @property
    def iam_instance_profile(self):
        return self.launch_spec["IamInstanceProfile"]

    def save(self):
        fw = open(self.filename, "w")
        s = json.dumps(self.spec, indent=4, sort_keys=True)
        # Clear the trailing spaces
        print("\n".join(x.rstrip() for x in s.splitlines()), file=fw)
        fw.close()

    def save_instance_id(self, instance_id, private_id_address):
        self.spec["InstanceId"] = instance_id
        self.spec["PrivateIpAddress"] = private_id_address
        self.save()

    def save_image_id(self, image_id):
        self.spec["LaunchSpec"]["ImageId"] = image_id
        self.save()


def main():

    actions = (
        ("cp", "copy files with support for wildcards"),
        ("ls", "list files with support for wildcards"),
        ("rm", "remove files with support for wildcards"),
        ("role", "change aws role"),
        ("start", "start ec2 instance"),
        ("stop", "stop ec2 instance"),
        ("ip", "describe current instance"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def ip(args):
    """
    %prog ip

    Show current IP address from JSON settings.
    """
    p = OptionParser(ip.__doc__)
    if len(args) != 0:
        sys.exit(not p.print_help())

    s = InstanceSkeleton()
    print("IP address:", s.private_ip_address, file=sys.stderr)
    print("Instance type:", s.instance_type, file=sys.stderr)


def start(args):
    """
    %prog start

    Launch ec2 instance through command line.
    """
    p = OptionParser(start.__doc__)
    p.add_option(
        "--ondemand",
        default=False,
        action="store_true",
        help="Do we want a more expensive on-demand instance",
    )
    p.add_option("--profile", default="mvrad-datasci-role", help="Profile name")
    p.add_option("--price", default=4.0, type=float, help="Spot price")
    opts, args = p.parse_args(args)

    if len(args) != 0:
        sys.exit(not p.print_help())

    role(["htang"])
    session = boto3.Session(profile_name=opts.profile)
    client = session.client("ec2")
    s = InstanceSkeleton()

    # Make sure the instance id is empty
    instance_id = s.instance_id
    if instance_id != "":
        logging.error("Instance exists {}".format(instance_id))
        sys.exit(1)

    launch_spec = s.launch_spec
    instance_id = ""

    if opts.ondemand:
        # Launch on-demand instance
        response = client.run_instances(
            BlockDeviceMappings=s.block_device_mappings,
            MaxCount=1,
            MinCount=1,
            ImageId=s.image_id,
            InstanceType=s.instance_type,
            KeyName=s.key_name,
            Placement={"AvailabilityZone": s.availability_zone},
            SecurityGroupIds=s.security_group_ids,
            SubnetId=s.subnet_id,
            EbsOptimized=s.ebs_optimized,
            IamInstanceProfile=s.iam_instance_profile,
        )
        instance_id = response["Instances"][0]["InstanceId"]

    else:
        # Launch spot instance
        response = client.request_spot_instances(
            SpotPrice=str(opts.price),
            InstanceCount=1,
            Type="one-time",
            AvailabilityZoneGroup=s.availability_zone,
            LaunchSpecification=launch_spec,
        )

        request_id = response["SpotInstanceRequests"][0]["SpotInstanceRequestId"]
        print("Request id {}".format(request_id), file=sys.stderr)

        while not instance_id:
            response = client.describe_spot_instance_requests(
                SpotInstanceRequestIds=[request_id]
            )
            if "InstanceId" in response["SpotInstanceRequests"][0]:
                instance_id = response["SpotInstanceRequests"][0]["InstanceId"]
            else:
                logging.debug("Waiting to be fulfilled ...")
                time.sleep(10)

    # Check if the instance is running
    print("Instance id {}".format(instance_id), file=sys.stderr)
    status = ""
    while status != "running":
        logging.debug("Waiting instance to run ...")
        time.sleep(3)
        response = client.describe_instance_status(InstanceIds=[instance_id])
        if len(response["InstanceStatuses"]) > 0:
            status = response["InstanceStatuses"][0]["InstanceState"]["Name"]

    # Tagging
    name = "htang-lx-ondemand" if opts.ondemand else "htang-lx-spot"
    response = client.create_tags(
        Resources=[instance_id],
        Tags=[
            {"Key": k, "Value": v}
            for k, v in {
                "Name": name,
                "owner": "htang",
                "project": "mv-bioinformatics",
            }.items()
        ],
    )

    # Attach working volumes
    volumes = s.volumes
    for volume in volumes:
        response = client.attach_volume(
            VolumeId=volume["VolumeId"], InstanceId=instance_id, Device=volume["Device"]
        )

    # Save instance id and ip
    response = client.describe_instances(InstanceIds=[instance_id])
    ip_address = response["Reservations"][0]["Instances"][0]["PrivateIpAddress"]
    print("IP address {}".format(ip_address), file=sys.stderr)

    s.save_instance_id(instance_id, ip_address)


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

    role(["htang"])
    session = boto3.Session(profile_name=opts.profile)
    client = session.client("ec2")
    s = InstanceSkeleton()

    # Make sure the instance id is NOT empty
    instance_id = s.instance_id
    if instance_id == "":
        logging.error("Cannot find instance_id {}".format(instance_id))
        sys.exit(1)

    block_device_mappings = []
    for volume in s.volumes:
        block_device_mappings.append({"DeviceName": volume["Device"], "NoDevice": ""})

    new_image_name = "htang-dev-{}-{}".format(timestamp(), int(time.time()))
    response = client.create_image(
        InstanceId=instance_id,
        Name=new_image_name,
        BlockDeviceMappings=block_device_mappings,
    )
    print(response, file=sys.stderr)
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
    old_snapshot_id = response["Images"][0]["BlockDeviceMappings"][0]["Ebs"][
        "SnapshotId"
    ]
    response = client.deregister_image(ImageId=old_image_id)
    print(response, file=sys.stderr)
    response = client.delete_snapshot(SnapshotId=old_snapshot_id)
    print(response, file=sys.stderr)
    response = client.terminate_instances(InstanceIds=[instance_id])
    print(response, file=sys.stderr)

    # Save new image id
    s.save_image_id(new_image_id)
    s.save_instance_id("", "")


def glob_s3(store, keys=None, recursive=False):
    store, cards = store.rsplit("/", 1)
    contents = ls_s3(store, recursive=recursive)
    if keys:
        filtered = [x for x in contents if op.basename(x).split(".")[0] in keys]
    else:
        filtered = fnmatch.filter(contents, cards)

    if recursive:
        store = "s3://" + store.replace("s3://", "").split("/")[0]

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

    (store,) = args
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
    p.add_option(
        "--force", default=False, action="store_true", help="Force overwrite if exists"
    )
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
        if isinstance(c, str):
            oc = op.basename(c)
            tc = op.join(folder, oc)
        else:
            if len(c) == 2:
                c, tc = c
            else:
                (c,) = c
                tc = op.basename(c)
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
    p.add_option(
        "--recursive", default=False, action="store_true", help="Recursive search"
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (store,) = args
    keys = opts.keys
    if keys:
        keys = SetFile(keys)
    print("\n".join(glob_s3(store, keys=keys, recursive=opts.recursive)))


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
    contents = []
    for row in popen(cmd):
        f = row.split()[-1]
        f = op.join(s3_store_obj_name, f)
        contents.append(f)

    if recursive:
        que = [x for x in contents if x.endswith("/")]
        while que:
            f = que.pop(0).rstrip("/")
            contents += ls_s3(f, recursive=True)

    return contents


def check_exists_s3(s3_store_obj_name):
    s3_store_obj_name = s3ify(s3_store_obj_name)
    cmd = "aws s3 ls {0} | wc -l".format(s3_store_obj_name)
    counts = int(popen(cmd).read())
    return counts != 0


def aws_configure(profile, key, value):
    sh("aws configure set profile.{0}.{1} {2}".format(profile, key, value))


def role(args):
    """
    %prog role htang

    Change aws role.
    """
    (
        src_acct,
        src_username,
        dst_acct,
        dst_role,
    ) = "205134639408 htang 114692162163 mvrad-datasci-role".split()

    p = OptionParser(role.__doc__)
    p.add_option("--profile", default="mvrad-datasci-role", help="Profile name")
    p.add_option(
        "--device",
        default="arn:aws:iam::" + src_acct + ":mfa/" + src_username,
        metavar="arn:aws:iam::123456788990:mfa/dudeman",
        help="The MFA Device ARN. This value can also be "
        "provided via the environment variable 'MFA_DEVICE' or"
        " the ~/.aws/credentials variable 'aws_mfa_device'.",
    )
    p.add_option(
        "--duration",
        type=int,
        default=3600,
        help="The duration, in seconds, that the temporary "
        "credentials should remain valid. Minimum value: "
        "900 (15 minutes). Maximum: 129600 (36 hours). "
        "Defaults to 43200 (12 hours), or 3600 (one "
        "hour) when using '--assume-role'. This value "
        "can also be provided via the environment "
        "variable 'MFA_STS_DURATION'. ",
    )
    p.add_option(
        "--assume-role",
        "--assume",
        default="arn:aws:iam::" + dst_acct + ":role/" + dst_role,
        metavar="arn:aws:iam::123456788990:role/RoleName",
        help="The ARN of the AWS IAM Role you would like to "
        "assume, if specified. This value can also be provided"
        " via the environment variable 'MFA_ASSUME_ROLE'",
    )
    p.add_option(
        "--role-session-name",
        help="Friendly session name required when using --assume-role",
        default=getpass.getuser(),
    )
    p.add_option(
        "--force",
        help="Refresh credentials even if currently valid.",
        action="store_true",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    # Use a config to check the expiration of session token
    config = get_config(AWS_CREDS_PATH)
    validate(opts, config)


def validate(args, config):
    """Validate if the config file is properly structured"""
    profile = args.profile
    if not args.profile:
        if os.environ.get("AWS_PROFILE"):
            args.profile = os.environ.get("AWS_PROFILE")
        else:
            args.profile = "default"

    if args.assume_role:
        role_msg = "with assumed role: %s" % (args.assume_role,)
    elif config.has_option(args.profile, "assumed_role_arn"):
        role_msg = "with assumed role: %s" % (
            config.get(args.profile, "assumed_role_arn")
        )
    else:
        role_msg = ""
    logging.info("Validating credentials for profile: %s %s" % (profile, role_msg))
    reup_message = "Obtaining credentials for a new role or profile."

    try:
        key_id = config.get(profile, "aws_access_key_id")
        access_key = config.get(profile, "aws_secret_access_key")
    except NoSectionError:
        log_error_and_exit(
            "Credentials session '[%s]' is missing. "
            "You must add this section to your credentials file "
            "along with your long term 'aws_access_key_id' and "
            "'aws_secret_access_key'" % (profile,)
        )
    except NoOptionError as e:
        log_error_and_exit(e)

    # get device from param, env var or config
    if not args.device:
        if os.environ.get("MFA_DEVICE"):
            args.device = os.environ.get("MFA_DEVICE")
        elif config.has_option(profile, "aws_mfa_device"):
            args.device = config.get(profile, "aws_mfa_device")
        else:
            log_error_and_exit(
                "You must provide --device or MFA_DEVICE or set "
                '"aws_mfa_device" in ".aws/credentials"'
            )

    # get assume_role from param or env var
    if not args.assume_role:
        if os.environ.get("MFA_ASSUME_ROLE"):
            args.assume_role = os.environ.get("MFA_ASSUME_ROLE")
        elif config.has_option(profile, "assume_role"):
            args.assume_role = config.get(profile, "assume_role")

    # get duration from param, env var or set default
    if not args.duration:
        if os.environ.get("MFA_STS_DURATION"):
            args.duration = int(os.environ.get("MFA_STS_DURATION"))
        else:
            args.duration = 3600 if args.assume_role else 43200

    # If this is False, only refresh credentials if expired. Otherwise
    # always refresh.
    force_refresh = False

    # Validate presence of profile-term section
    if not config.has_section(profile):
        config.add_section(profile)
        force_refresh = True
    # Validate option integrity of profile section
    else:
        required_options = [
            "assumed_role",
            "aws_access_key_id",
            "aws_secret_access_key",
            "aws_session_token",
            "aws_security_token",
            "expiration",
        ]
        try:
            short_term = {}
            for option in required_options:
                short_term[option] = config.get(profile, option)
        except NoOptionError:
            logging.warning(
                "Your existing credentials are missing or invalid, "
                "obtaining new credentials."
            )
            force_refresh = True

        try:
            current_role = config.get(profile, "assumed_role_arn")
        except NoOptionError:
            current_role = None

        if args.force:
            logging.info("Forcing refresh of credentials.")
            force_refresh = True
        # There are not credentials for an assumed role,
        # but the user is trying to assume one
        elif current_role is None and args.assume_role:
            logging.info(reup_message)
            force_refresh = True
        # There are current credentials for a role and
        # the role arn being provided is the same.
        elif (
            current_role is not None
            and args.assume_role
            and current_role == args.assume_role
        ):
            pass
        # There are credentials for a current role and the role
        # that is attempting to be assumed is different
        elif (
            current_role is not None
            and args.assume_role
            and current_role != args.assume_role
        ):
            logging.info(reup_message)
            force_refresh = True
        # There are credentials for a current role and no role arn is
        # being supplied
        elif current_role is not None and args.assume_role is None:
            logging.info(reup_message)
            force_refresh = True

    should_refresh = True

    # Unless we're forcing a refresh, check expiration.
    if not force_refresh:
        exp = datetime.strptime(config.get(profile, "expiration"), "%Y-%m-%d %H:%M:%S")
        diff = exp - datetime.utcnow()
        if diff.total_seconds() <= 0:
            logging.info("Your credentials have expired, renewing.")
        else:
            should_refresh = False
            logging.info(
                "Your credentials are still valid for %s seconds"
                " they will expire at %s" % (diff.total_seconds(), exp)
            )

    if should_refresh:
        get_credentials(profile, args, config)


def get_credentials(profile, args, config):
    mfa_token = console.input(
        "Enter AWS MFA code for device [%s] "
        "(renewing for %s seconds): " % (args.device, args.duration)
    )

    boto3.setup_default_session(profile_name="default")
    client = boto3.client("sts")

    if args.assume_role:

        logging.info(
            "Assuming Role - Profile: %s, Role: %s, Duration: %s",
            profile,
            args.assume_role,
            args.duration,
        )

        try:
            print((args.assume_role, args.role_session_name, args.device, mfa_token))
            response = client.assume_role(
                RoleArn=args.assume_role,
                RoleSessionName=args.role_session_name,
                SerialNumber=args.device,
                TokenCode=mfa_token,
            )
        except ClientError as e:
            log_error_and_exit(
                "An error occured while calling assume role: {}".format(e)
            )
        except ParamValidationError:
            log_error_and_exit("Token must be six digits")

        config.set(
            profile,
            "assumed_role",
            "True",
        )
        config.set(
            profile,
            "assumed_role_arn",
            args.assume_role,
        )
    else:
        logging.info(
            "Fetching Credentials - Profile: %s, Duration: %s", profile, args.duration
        )
        try:
            response = client.get_session_token(
                DurationSeconds=args.duration,
                SerialNumber=args.device,
                TokenCode=mfa_token,
            )
        except ClientError as e:
            log_error_and_exit(
                "An error occured while calling assume role: {}".format(e)
            )
        except ParamValidationError:
            log_error_and_exit("Token must be six digits")

        config.set(
            profile,
            "assumed_role",
            "False",
        )
        config.remove_option(profile, "assumed_role_arn")

    # aws_session_token and aws_security_token are both added
    # to support boto and boto3
    options = [
        ("aws_access_key_id", "AccessKeyId"),
        ("aws_secret_access_key", "SecretAccessKey"),
        ("aws_session_token", "SessionToken"),
        ("aws_security_token", "SessionToken"),
    ]

    for option, value in options:
        config.set(profile, option, response["Credentials"][value])
    # Save expiration individiually, so it can be manipulated
    config.set(
        profile,
        "expiration",
        response["Credentials"]["Expiration"].strftime("%Y-%m-%d %H:%M:%S"),
    )
    with open(AWS_CREDS_PATH, "w") as configfile:
        config.write(configfile)

    logging.info(
        "Success! Your credentials will expire in %s seconds at: %s"
        % (args.duration, response["Credentials"]["Expiration"])
    )


def log_error_and_exit(message):
    """Log an error message and exit with error"""
    logging.error(message)
    sys.exit(1)


if __name__ == "__main__":
    main()
