#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
AWS-related methods.
"""

import os.path as op
import sys
import fnmatch

import boto3
from multiprocessing import Pool

from jcvi.formats.base import SetFile
from jcvi.apps.base import OptionParser, ActionDispatcher, popen, sh


def main():

    actions = (
        ('cp', 'copy files with support for wildcards'),
        ('ls', 'list files with support for wildcards'),
        ('rm', 'remove files with support for wildcards'),
        ('role', 'change aws role'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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
        if len(c) == 1:
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
