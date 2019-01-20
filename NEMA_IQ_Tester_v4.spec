# -*- mode: python -*-

block_cipher = None


a = Analysis(['NEMA_IQ_Tester_v4.py'],
             pathex=['E:\\NEMA_IQ'],
             binaries=None,
             datas=None,
             hiddenimports=['scipy.integrate', 'scipy.integrate.quadpack', 'scipy.integrate._vode'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='NEMA_IQ_Tester_v4',
          debug=False,
          strip=False,
          upx=True,
          console=True )
