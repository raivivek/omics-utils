#! /usr/bin/env python3


from argh import *

@arg('-a', 'analysis', help='Where the ANALYSIS goes')
@arg('-c', 'control', help='Where project CONTROL is setup')
def init(control, analysis):
    """ Setup a skeleton project at the specified ANALYSIS directory with
    CONTROL setup at a different (Version Controlled) location. """

    raise NotImplementedError


def save(target):
    """ Commit all changes and push to Git remote """

    raisee NotImplementedError


def snakemake():


if __name__ == '__main__':
    parser = ArghParser()
    parser.add_commands([init])
    parser.dispatch()
