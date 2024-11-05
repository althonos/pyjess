set(JESS_SOURCES
    Annulus.c
    Annulus.h
    Atom.c
    Atom.h
    Jess.c
    Jess.h
    Join.c
    Join.h
    KdTree.c
    KdTree.h
    Molecule.c
    Molecule.h
    Region.c
    Region.h
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
target_compile_definitions(jess PUBLIC -D_GNU_SOURCE)
set_property(TARGET jess PROPERTY C_STANDARD 99)