REM This requires mingw-w64 installed with x86_64, win32, and seh (though seh doesn't matter).

"C:\Program Files\mingw-w64\x86_64-8.1.0-win32-seh-rt_v6-rev0\mingw64\bin\gcc.exe" -shared -Wl,-soname,engine -o engine.dll -fPIC engine.c

@pause
