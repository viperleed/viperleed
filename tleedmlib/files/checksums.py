# -*- coding: utf-8 -*-
"""
@author: Alexander M. Imre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! MAKE SURE THIS FILE HAS NOT BEEN TAMPERD WITH !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Contains checksums for TensErLEED files.
Must be updated for each new version of TensErLEED. Use function 
_generate_checksums_for_dir() from TL_base.py to quickly generate 
the required checksums.
"""

# permissible checksums for various TensErLEED version
# implemented as a dict of dicts - outer key is TL version, inner key is filename
# allowed checksums for each version are given in a tuple

_TL_1_6_INPUT_FILES = {
    'TensErLEED-v1.6/src/superpos.v1.6.f':
        ('03a8546e3929f3ee07aea124f0d5fb8393332a0950c99cca830383d68c45f765', ),
    'TensErLEED-v1.6/src/search.mpi.f':
        ('95f6fd52e60cbfe25f1f095a77f91858d530c73301ead17a5d3b9c1e11435233', ),
    'TensErLEED-v1.6/src/ref-calc.v1.6.f':
        ('e073f4230f2540b070e2989c6f1e52ad524fe13e1774b90bea305689542593bb', ),
    'TensErLEED-v1.6/src/delta.v1.6.f':
        ('bc9a5524f19b1cb62d44264158a253a7e2719af51e39c9cde5eb7c81cdf3d7a1', ),
    'TensErLEED-v1.6/src/rfactor.v1.6.f':
        ('49d87897a07c55d9782d85860cd279b7dd3d1592ba80ae278a67908eb60d8549', ),
    'TensErLEED-v1.6/src/search.v1.6.f':
        ('d68b53a47beeb16fb002a9bbacc347fe7685175e91e5393533c9a76dd39ed782', ),
    'TensErLEED-v1.6/src/rfactorforerror.v1.6.f':
        ('60170fae9632f17b3c7454029e5cfdcf2815fe9e0cd209707774d76d0919f581', ),
    'TensErLEED-v1.6/lib/lib.delta.v1.6.f':
        ('f3098275842a221ec865f03314f1f6012f3209900714fb55142c07ec96c4cdb6', ),
    'TensErLEED-v1.6/lib/lib.search.v1.6.f':
        ('8ac98c917042feca21e81854ea0da1041768e33ea4402a4969af89c3da55c2d5', ),
    'TensErLEED-v1.6/lib/lib.search.mpi.f':
        ('f5a8ac72ad4ec220d0be12912f6cd6f91d92f6da539ded65597db9433902b555', ),
    'TensErLEED-v1.6/lib/rfacsb.v1.6.f':
        ('66f4950480dfc7f4f3cfb99c310217e27559b2a37992365752ae9bc99eaa218c', ),
    'TensErLEED-v1.6/lib/lib.superpos.v1.6.f':
        ('1548ab1c020fffc912b7c92b2bda05d040cad60ddf710fc7350374a554f1d6b7', ),
    'TensErLEED-v1.6/lib/lib.tleed.v1.6.f':
        ('7e46e26f12dc7f4b82c227b865b08ef1962c58ed4d6ff3d658fce7f0827b8eec', ),
}

_TL_1_61_INPUT_FILES = {
    'TensErLEED-v1.61/src/superpos.v1.6.f':
        ('03a8546e3929f3ee07aea124f0d5fb8393332a0950c99cca830383d68c45f765', ),
    'TensErLEED-v1.61/src/ref-calc.v1.6.f':
        ('e073f4230f2540b070e2989c6f1e52ad524fe13e1774b90bea305689542593bb', ),
    'TensErLEED-v1.61/src/rfactor.v1.6.f':
        ('49d87897a07c55d9782d85860cd279b7dd3d1592ba80ae278a67908eb60d8549', ),
    'TensErLEED-v1.61/src/search.v1.6.f':
        ('d68b53a47beeb16fb002a9bbacc347fe7685175e91e5393533c9a76dd39ed782', ),
    'TensErLEED-v1.61/src/delta.v1.61.f':
        ('1acdaaf0bd5984b0b477d11f89ff00c884fd231a26a0d2aa2e7909c2c3c25877', ),
    'TensErLEED-v1.61/src/rfactorforerror.v1.6.f':
        ('60170fae9632f17b3c7454029e5cfdcf2815fe9e0cd209707774d76d0919f581', ),
    'TensErLEED-v1.61/src/search.mpi.f90':
        ('55ffec2df2a9d573ccdfb6b1fa5e7c9706158aa96ac3618e68dd277a37d22a7c', ),
    'TensErLEED-v1.61/lib/lib.search.v1.6.f':
        ('8ac98c917042feca21e81854ea0da1041768e33ea4402a4969af89c3da55c2d5', ),
    'TensErLEED-v1.61/lib/lib.search.mpi.f':
        ('ae9bf35bc5270bf50ddad66ce3465fa797324ffaa323c0b77fa7c4bfa65cdbec', ),
    'TensErLEED-v1.61/lib/rfacsb.v1.6.f':
        ('66f4950480dfc7f4f3cfb99c310217e27559b2a37992365752ae9bc99eaa218c', ),
    'TensErLEED-v1.61/lib/lib.superpos.v1.6.f':
        ('1548ab1c020fffc912b7c92b2bda05d040cad60ddf710fc7350374a554f1d6b7', ),
    'TensErLEED-v1.61/lib/lib.delta.v1.61.f':
        ('33f96300f6e823f796668b89ab9e8821460efb7b716381e76a6a1be15be48a63', ),
    'TensErLEED-v1.61/lib/intarr_hashing.f90':
        ('0eafda23bcd4378daf106016c7e35334fad49bb582b0c57b226ac8449c98173c', ),
    'TensErLEED-v1.61/lib/lib.tleed.v1.6.f':
        ('7e46e26f12dc7f4b82c227b865b08ef1962c58ed4d6ff3d658fce7f0827b8eec', ),
}

_TL_1_71_INPUT_FILES = {
    'TensErLEED-v1.71/src/superpos.v1.71.f':
        ('bc3211d313448073d7a0d16f289592479b6d99994facb16cdb2d98b3ad143e24', ),
    'TensErLEED-v1.71/src/rfactor.v1.71.f':
        ('b770c1e5bf0e28bd6c6c55d19c71a66bd3a64307872e70e8702592c6d9463c94', ),
    'TensErLEED-v1.71/src/search.mpi.v1.71.f90':
        ('2c4e00cdb9eb0beb5705a6fc04d871221969120bebd6e548c72b84f8068c0919', ),
    'TensErLEED-v1.71/src/delta.v1.71.f':
        ('eadbfb0595524ee88de5f241b040b082696615fdf8c8e2b842f0595378ebebbe', ),
    'TensErLEED-v1.71/src/ref-calc.v1.71.f':
        ('7a83754d75a98733fd3679dc62d587870b52282a0f69f07baf4f24cb1ebd07a6', ),
    'TensErLEED-v1.71/src/search.v1.71.f90':
        ('bb9ba7c6a959766a16db7b4b3566265542f4e2c95afde52f636c4a1351ab0e4c', ),
    'TensErLEED-v1.71/lib/lib.superpos.v1.71.f':
        ('e685755c7176233e01cd78ed7c6790d648007664d14e76d8657c774614fd64af', ),
    'TensErLEED-v1.71/lib/lib.search.v1.71.f':
        ('8c9634e34bca746817a6d2bf8ce324c9fe8bba63fd27eb2df64245cd98e44f1d', ),
    'TensErLEED-v1.71/lib/lib.tleed.v1.71.f':
        ('12a496e924a2a3c09d176353f15df6ad8da0f701635d122f2ec7af6010bf0107', ),
    'TensErLEED-v1.71/lib/rfacsb.v1.71.f':
        ('d6d5fdaa5c96a84e35514b9bf246a173244d8a661dd72da2909d9c225d1e237b', ),
    'TensErLEED-v1.71/lib/intarr_hashing.f90':
        ('0eafda23bcd4378daf106016c7e35334fad49bb582b0c57b226ac8449c98173c', ),
    'TensErLEED-v1.71/lib/lib.delta.v1.71.f':
        ('2ee77a2fdf4abdb90e55a6107d8948576d60f91de4eacf3c0bc7fc53a8702085', ),
    'TensErLEED-v1.71/lib/lib.search.mpi.v1.71.f':
        ('f761a8e66337ea809951051dae5128637ce8135d55dba8b8ce88970b5c700067', ),
}

_TL_1_72_INPUT_FILES = {
    'TensErLEED-v1.72/src/rfactor.v1.72.f':
        ('1e066ee91474b9556e79aefa844670833bc4559b9799b4c8d43e7b94e3ca455f', ),
    'TensErLEED-v1.72/src/superpos.v1.72.f':
        ('bc3211d313448073d7a0d16f289592479b6d99994facb16cdb2d98b3ad143e24', ),
    'TensErLEED-v1.72/src/search.mpi.v1.72.f90':
        ('2c4e00cdb9eb0beb5705a6fc04d871221969120bebd6e548c72b84f8068c0919', ),
    'TensErLEED-v1.72/src/search.v1.72.f90':
        ('bb9ba7c6a959766a16db7b4b3566265542f4e2c95afde52f636c4a1351ab0e4c', ),
    'TensErLEED-v1.72/src/delta.v1.72.f':
        ('eadbfb0595524ee88de5f241b040b082696615fdf8c8e2b842f0595378ebebbe', ),
    'TensErLEED-v1.72/src/ref-calc.v1.72.f':
        ('bf5eb2ed41ec690ffd3203a2ce8946a5eb33825622a25ade6963e2b9ba96b0e1', ),
    'TensErLEED-v1.72/lib/lib.superpos.v1.72.f':
        ('e685755c7176233e01cd78ed7c6790d648007664d14e76d8657c774614fd64af', ),
    'TensErLEED-v1.72/lib/lib.search.v1.72.f':
        ('15760921e03c562cfa01084745b0ab73950a2e988b0e42f13f1c09ef13ba8d9e', ),
    'TensErLEED-v1.72/lib/rfacsb.v1.72.f':
        ('978fa671cbfdd8905cb9e799bc30e2d79de91398fd1b7ccbb13e794aa6cf0106', ),
    'TensErLEED-v1.72/lib/lib.delta.v1.72.f':
        ('2ee77a2fdf4abdb90e55a6107d8948576d60f91de4eacf3c0bc7fc53a8702085', ),
    'TensErLEED-v1.72/lib/lib.search.mpi.v1.72.f':
        ('f761a8e66337ea809951051dae5128637ce8135d55dba8b8ce88970b5c700067', ),
    'TensErLEED-v1.72/lib/intarr_hashing.f90':
        ('0eafda23bcd4378daf106016c7e35334fad49bb582b0c57b226ac8449c98173c', ),
    'TensErLEED-v1.72/lib/lib.tleed.v1.72.f':
        ('d76968ec43b7697cd517a65856cfc4cadcce2f09f0d6ecc57700cc70f2617520', ),
}

_TL_1_73_INPUT_FILES = {
    'TensErLEED-v1.73/src/ref-calc.f':
        ('297597429e66505d0ece90bcc4c0b9da82745f28b96c358066981bfbfae27237', ),
    'TensErLEED-v1.73/src/delta.f':
        ('eadbfb0595524ee88de5f241b040b082696615fdf8c8e2b842f0595378ebebbe', ),
    'TensErLEED-v1.73/src/search.mpi.f90':
        ('1827fb5ae553e52becf4365a0cceb8a6acbdf96b97692662d2d702941862ef3a', ),
    'TensErLEED-v1.73/src/rfactor.f':
        ('62bfe5bf8543d47b2e460b04244e0593ef98524746c22e4646123c4e6a172541', ),
    'TensErLEED-v1.73/src/superpos.f':
        ('bc3211d313448073d7a0d16f289592479b6d99994facb16cdb2d98b3ad143e24', ),
    'TensErLEED-v1.73/src/search.f90':
        ('bb9ba7c6a959766a16db7b4b3566265542f4e2c95afde52f636c4a1351ab0e4c', ),
    'TensErLEED-v1.73/lib/lib.search.f':
        ('6da7874dda85168d972c6e04cc6c134ea5819b67b511510824d306854641cb3d', ),
    'TensErLEED-v1.73/lib/lib.delta.f':
        ('2ee77a2fdf4abdb90e55a6107d8948576d60f91de4eacf3c0bc7fc53a8702085', ),
    'TensErLEED-v1.73/lib/lib.search.mpi.f':
        ('e34b0c8935c5723b9df88d723bbaa2d8f2bde6c2b6ab2a0bbf5d0aed1fce2fd3', ),
    'TensErLEED-v1.73/lib/lib.tleed.f':
        ('fc800ff5f9d97c633700736675933d784a52291890aebed5f57e0011425d63c3', ),
    'TensErLEED-v1.73/lib/rfacsb.f':
        ('978fa671cbfdd8905cb9e799bc30e2d79de91398fd1b7ccbb13e794aa6cf0106', ),
    'TensErLEED-v1.73/lib/lib.superpos.f':
        ('e685755c7176233e01cd78ed7c6790d648007664d14e76d8657c774614fd64af', ),
    'TensErLEED-v1.73/lib/intarr_hashing.f90':
        ('0eafda23bcd4378daf106016c7e35334fad49bb582b0c57b226ac8449c98173c', ),
}

INPUT_FILE_VERSIONS = {
    '1.6': _TL_1_6_INPUT_FILES,
    '1.61': _TL_1_61_INPUT_FILES,
    '1.71': _TL_1_71_INPUT_FILES,
    '1.72': _TL_1_72_INPUT_FILES,
    '1.73': _TL_1_73_INPUT_FILES,
}
