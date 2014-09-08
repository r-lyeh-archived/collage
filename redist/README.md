collage/redist
==============

- This optional folder is used to regenerate the amalgamated distribution. Do not include it into your project.
- Regenerate the distribution by typing the following lines:
```
move /y collage.hpp ..
deps\Amalgamate.exe -w "*.*pp;*.c;*.h" collage.cpp ..\collage.cpp
deps\fart.exe -- ..\collage.cpp "#line" "//#line"
deps\fart.exe -- ..\collage.cpp "#pragma once" "//#pragma once"
copy ..\collage.hpp /y
```
