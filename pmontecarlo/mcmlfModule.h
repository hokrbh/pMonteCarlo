#ifndef MCMLFMODULE_H
#define MCMLFMODULE_H

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

int parseInputString(const char* cfgString, PROP *cfg);

#endif
