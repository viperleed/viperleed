# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['E:\\Work\\Vienna\\SPECSLab\\Technical\\Upgrades\\Python LEED-IV\\000-git\\viperleed\\gui.py'],
             pathex=['E:\\Work\\Vienna\\SPECSLab\\Technical\\Upgrades\\Python LEED-IV\\000-git\\00-GUI_builds\\v0.1.0\\'],  # perhaps need another \\ at end?
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
          name='viperleed_gui_v0.1.0',  #TODO: get version automatically
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True )
