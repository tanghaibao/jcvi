#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
AWS-related methods.
"""

import os.path as op
import sys

import boto3

from jcvi.apps.base import OptionParser, ActionDispatcher, sh


def main():

    actions = (
        ('role', 'change aws role'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def push_to_s3(s3_store, dir_name):
    cmd = "aws s3 sync %s s3://%s/%s --sse"%(dir_name, s3_store, dir_name)
    sh(cmd)


def pull_from_s3(s3_store):
    file_name = s3_store.split("/")[-1]
    cmd = "aws s3 cp s3://%s %s --sse"%(s3_store, file_name)
    sh(cmd)
    return op.abspath(file_name)


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
