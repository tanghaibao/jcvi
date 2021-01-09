#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# ez_setup.py
# utils
#
# Created by Haibao Tang on 11/24/20
# Copyright © 2021 Haibao Tang. All rights reserved.
#

"""
Identify the best downloading mechanism for a given URL.
Credits: https://pypi.org/project/ez_setup/
"""


import os
import platform
import subprocess

from urllib.request import urlopen


def download_file_powershell(url, target, cookies=None):
    """
    Download the file at url to target using Powershell (which will validate
    trust). Raise an exception if the command cannot complete.
    """
    if cookies:
        raise NotImplementedError
    target = os.path.abspath(target)
    cmd = [
        "powershell",
        "-Command",
        f"(new-object System.Net.WebClient).DownloadFile({url}, {target})",
    ]
    subprocess.check_call(cmd)


def has_powershell():
    if platform.system() != "Windows":
        return False
    cmd = ["powershell", "-Command", "echo test"]
    devnull = open(os.path.devnull, "wb")
    try:
        try:
            subprocess.check_call(cmd, stdout=devnull, stderr=devnull)
        except FileNotFoundError:
            return False
    finally:
        devnull.close()
    return True


download_file_powershell.viable = has_powershell


def download_file_curl(url, target, cookies=None):
    cmd = ["curl", url, "--output", target]
    # https://github.com/tanghaibao/jcvi/issues/307
    # When downloading Phytozome directory listing, there are multiple redirects
    # before we hit the index page. Natually we'd follow the redirects, similar
    # to the default behavior of wget
    cmd += ["-L"]  # follow redirect
    if url.startswith("ftp:"):
        cmd += ["-P", "-"]
    if cookies:
        cmd += ["-b", cookies]
    subprocess.check_call(cmd)


def has_curl():
    cmd = ["curl", "--version"]
    devnull = open(os.path.devnull, "wb")
    try:
        try:
            subprocess.check_call(cmd, stdout=devnull, stderr=devnull)
        except FileNotFoundError:
            return False
    finally:
        devnull.close()
    return True


download_file_curl.viable = has_curl


def download_file_wget(url, target, cookies=None):
    cmd = ["wget", url, "--output-document", target]
    cmd += ["--no-check-certificate"]
    if url.startswith("ftp:"):
        cmd += ["--passive-ftp"]
    if cookies:
        cmd += ["--load-cookies", cookies]
    subprocess.check_call(cmd)


def has_wget():
    cmd = ["wget", "--version"]
    devnull = open(os.path.devnull, "wb")
    try:
        try:
            subprocess.check_call(cmd, stdout=devnull, stderr=devnull)
        except FileNotFoundError:
            return False
    finally:
        devnull.close()
    return True


download_file_wget.viable = has_wget


def download_file_insecure(url, target, cookies=None):
    """
    Use Python to download the file, even though it cannot authenticate the
    connection.
    """
    if cookies:
        raise NotImplementedError
    src = dst = None
    try:
        src = urlopen(url)
        # Read/write all in one block, so we don't create a corrupt file
        # if the download is interrupted.
        data = src.read()
        dst = open(target, "wb")
        dst.write(data)
    finally:
        if src:
            src.close()
        if dst:
            dst.close()


download_file_insecure.viable = lambda: True

ALL_DOWNLOADERS = [
    ("wget", download_file_wget),
    ("curl", download_file_curl),
    ("powershell", download_file_powershell),
    ("insecure", download_file_insecure),
]


def get_best_downloader(downloader=None):
    """Choose among a set of 4 popular downloaders, in the following order:
    - wget
    - curl
    - powershell
    - insecure (Python)

    Args:
        downloader (str, optional): Use a given downloader. One of wget|curl|powershell|insecure.
        Defaults to None.

    Returns:
        Download function: The downloader function that accepts as parameters url, target
        and cookies.
    """
    for dl_name, dl in ALL_DOWNLOADERS:
        if downloader and dl_name != downloader:
            continue
        if dl.viable():
            return dl
