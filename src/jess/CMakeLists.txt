# --- Check if GNU extensions are available ------------------------------------

include(CheckSymbolExists)

check_symbol_exists(alloca "alloca.h" HAVE_ALLOCA)
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    check_symbol_exists(__GNU_LIBRARY__ "stdio.h" GLIBC)
endif()

# --- Detect reentrant qsort ---------------------------------------------------

include(CheckSourceRuns)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    check_source_runs(C
        "
        #include <stdlib.h>

        static int
        cmpstringp(void* data, const void *p1, const void *p2)
        { return 0; }

        int
        main(int argc, const char** argv)
        { qsort_r(&argv[1], argc - 1, sizeof(char *), NULL, cmpstringp); }
        "
        HAVE_APPLE_QSORT_R
    )

elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    check_source_runs(C
        "
        #include <stdlib.h>

        static int
        cmpstringp(void* data, const void *p1, const void *p2)
        { return 0; }

        int
        main(int argc, const char** argv)
        { qsort_s(&argv[1], argc - 1, sizeof(char *), cmpstringp, NULL); }
        "
        HAVE_WIN32_QSORT_S
    )
else()
    set(SAFE_CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS}")
    if(GLIBC)
        set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS}" -D_GNU_SOURCE)
    endif()
    check_source_runs(C
        "
        #include <stdlib.h>

        static int
        cmpstringp(const void *p1, const void *p2, void* data)
        { return 0; }

        int
        main(int argc, const char** argv)
        { qsort_r(&argv[1], argc - 1, sizeof(char *), cmpstringp, NULL); }
        "
        HAVE_GNU_QSORT_R
    )
    set(CMAKE_REQUIRED_DEFINITIONS "${SAFE_CMAKE_REQUIRED_DEFINITIONS}")
endif()

check_source_runs(C
    "
    #include <stdlib.h>
    static __thread int i;
    int main(int argc, const char** argv) { i = 0; return i; }
    "
    HAVE_THREADLOCALSTORAGE
)

# --- Compile Jess -------------------------------------------------------------

set(JESS_SOURCES
    Annulus.c
    Annulus.h
    Atom.c
    Atom.h
    Box.c
    Box.h
    CandidateSet.c
    CandidateSet.h
    Jess.c
    Jess.h
    Join.c
    Join.h
    KdTree.c
    KdTree.h
    Molecule.c
    Molecule.h
    qselect.h
    Region.c
    Region.h
    ResIndex.c
    ResIndex.h
    Scanner.c
    Scanner.h
    Super.c
    Super.h
    Template.h
    TessAtom.c
    TessAtom.h
    TessTemplate.c
    TessTemplate.h
)

foreach(_file IN ITEMS ${JESS_SOURCES})

    if(EXISTS ${PROJECT_SOURCE_DIR}/patches/${_file}.patch)


        if(NOT EXISTS ${PROJECT_SOURCE_DIR}/vendor/jess/src/${_file})
            add_custom_command(
                OUTPUT
                    ${_file}
                COMMENT
                    "Creating ${_file}"
                COMMAND
                    cmake -E touch ${CMAKE_CURRENT_BINARY_DIR}/empty.txt
                COMMAND
                    ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                        --input ${CMAKE_CURRENT_BINARY_DIR}/empty.txt
                        --patch ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
                        --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
                DEPENDS
                    ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
            )
        else()
            add_custom_command(
                OUTPUT
                    ${_file}
                COMMENT
                    "Patching ${_file}"
                COMMAND
                    ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                        --input ${PROJECT_SOURCE_DIR}/vendor/jess/src/${_file}
                        --patch ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
                        --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
                DEPENDS
                    ${PROJECT_SOURCE_DIR}/vendor/jess/src/${_file}
                    ${PROJECT_SOURCE_DIR}/patches/${_file}.patch
            )
        endif()
    else()
        add_custom_command(
            OUTPUT
                ${_file}
            COMMAND
                cmake -E copy ${PROJECT_SOURCE_DIR}/vendor/jess/src/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/jess/src/${_file}
        )
    endif()
    set(JESS_PATCHED_SOURCES ${JESS_PATCHED_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endforeach()

add_library(jess STATIC ${JESS_PATCHED_SOURCES})
target_include_directories(jess PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
set_property(TARGET jess PROPERTY C_STANDARD 11)
if(GLIBC)
    target_compile_definitions(jess PRIVATE -D_GNU_SOURCE)
endif()
if(HAVE_ALLOCA)
    target_compile_definitions(jess PRIVATE -DHAVE_ALLOCA)
endif()
if(HAVE_APPLE_QSORT_R)
    target_compile_definitions(jess PRIVATE -DHAVE_APPLE_QSORT_R)
elseif(HAVE_GNU_QSORT_R)
    target_compile_definitions(jess PRIVATE -DHAVE_GNU_QSORT_R)
elseif(HAVE_WIN32_QSORT_S)
    target_compile_definitions(jess PRIVATE -DHAVE_WIN32_QSORT_S)
endif()
if(HAVE_THREADLOCALSTORAGE)
    target_compile_definitions(jess PRIVATE -DHAVE_THREADLOCALSTORAGE)
endif()
