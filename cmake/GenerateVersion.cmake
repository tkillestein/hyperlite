# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
# Generates versionmod.f90 from versionmod.f90.in at *build* time so the
# stamped git revision tracks HEAD without a reconfigure.  Invoked as
#   cmake -DSRC=... -DDST=... -DHYPERLITE_VER=... -DSOURCE_DIR=... -P GenerateVersion.cmake
# copy_if_different keeps the mtime stable when nothing changed, so an
# unchanged revision does not trigger a rebuild cascade.

execute_process(
  COMMAND git -C "${SOURCE_DIR}" describe --tags --always --dirty
  OUTPUT_VARIABLE HYPERLITE_GIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
  RESULT_VARIABLE _git_rc)
if(NOT _git_rc EQUAL 0 OR HYPERLITE_GIT STREQUAL "")
  set(HYPERLITE_GIT "unknown")
endif()

configure_file("${SRC}" "${DST}.tmp" @ONLY)
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different "${DST}.tmp" "${DST}")
file(REMOVE "${DST}.tmp")
