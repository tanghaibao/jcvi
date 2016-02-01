#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
AWS-related methods.
"""

import os.path as op
import sys

import boto3

from jcvi.apps.base import OptionParser, ActionDispatcher, popen, sh


def main():

    actions = (
        ('role', 'change aws role'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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


def pull_from_s3(s3_store, file_name=None):
    file_name = file_name or s3_store.split("/")[-1]
    s3_store = s3ify(s3_store)
    cmd = "aws s3 cp {0} {1} --sse".format(s3_store, file_name)
    sh(cmd)
    return op.abspath(file_name)


def ls_s3(s3_store_obj_name):
    s3_store_obj_name = s3ify(s3_store_obj_name)
    cmd = "aws s3 ls {0}/".format(s3_store_obj_name)
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
