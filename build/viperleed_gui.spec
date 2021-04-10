# -*- mode: python ; coding: utf-8 -*-

# Notes for running:
#
# py -3 -m PyInstaller viperleed_gui.spec --onefile
#
# IF IT GETS STUCK at PKG-00:
# remove the dist and build directories, and run from an admin cmd
# It can also help to rather run
#   %path_to_python_version%/python.exe -m PyInstaller viperleed_gui.spec --onefile

# Matplotlib is a mess.
# v3.4 renamed the _get_data_path() to get_data_path, so one needs
# to edit the PyInstaller/hooks/hook-matplotlib.py to replace
# the invalid name. Otherwise, the exe builds, but will not find the
# mpl-data directory in the extracted archive

block_cipher = None


a = Analysis(['E:\\Work\\Vienna\\SPECSLab\\Technical\\Upgrades\\Python LEED-IV\\000-git\\viperleed\\gui.py'],
             pathex=['E:\\Work\\Vienna\\SPECSLab\\Technical\\Upgrades\\Python LEED-IV\\000-git\\00-GUI_builds\\v0.3.0\\'],
             binaries=[],
             datas=[],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)

font_path = ("E:\\Work\\Vienna\\SPECSLab\\Technical\\Upgrades\\Python LEED-IV"
			 "\\000-git\\viperleed\\guilib\\fonts\\\\DejaVuSans.ttf")

font_files = [("guilib/fonts/DejaVuSans.ttf", font_path, "DATA"),
			  ("guilib/fonts/cmunrm.otf", font_path, "DATA")]

pyz = PYZ(a.pure, a.zipped_data,
          cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas + font_files,
          [],
          name='viperleed_gui_v0.3.0',  #TODO: get version automatically
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )
