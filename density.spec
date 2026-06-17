# -*- mode: python ; coding: utf-8 -*-

from PyInstaller.utils.hooks import collect_data_files, collect_submodules


thermo_datas = collect_data_files('thermo')
chemicals_datas = collect_data_files('chemicals')
thermo_hiddenimports = collect_submodules('thermo')
chemicals_hiddenimports = collect_submodules('chemicals')

a = Analysis(
    ['density.py'],
    pathex=[],
    binaries=[],
    datas=[
        ('lhv_data.db', '.'),
        ('gui.ui', '.'),
        *thermo_datas,
        *chemicals_datas,
    ],
    hiddenimports=thermo_hiddenimports + chemicals_hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='density ver3',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
