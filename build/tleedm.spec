# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['../tleedm.py'],
             datas=[('../tleedmlib/_checksums.dat', 'viperleed/tleedmlib/'), ('../tleedm.mplstyle', 'viperleed/')],
             hiddenimports=['libgcc','scipy.spacial._kdtree','scipy.spacial._cdktree','sklearn.utils._cython_blas', 'scipy.spatial.transform._rotation_groups', 'scipy.special.cython_special', 'sklearn.utils._weight_vector', 'sklearn.utils._typedefs', 'sklearn.utils._heap', 'sklearn.utils._sorting', 'sklearn.utils._vector_sentinel'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='tleedm',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )
