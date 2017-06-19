cmake_policy(SET CMP0020 NEW)

aux_source_directory(. SOURCE_FILE_LIST)
file (GLOB_RECURSE HEADER_FILES RELATIVE ${PROJECT_SOURCE_DIR} *.h *.hpp *.hxx)
file (GLOB_RECURSE RESOURCES RELATIVE ${PROJECT_SOURCE_DIR} *.qrc)
file (GLOB_RECURSE WINRC_FILES RELATIVE ${PROJECT_SOURCE_DIR} *.rc)
qt5_wrap_cpp (MOC_SOURCES ${HEADER_FILES})
qt5_add_resources (QRC_SOURCES ${RESOURCES})
qt5_create_translation (QM_FILE ${SOURCE_FILE_LIST} ${HEADER_FILES} ${TS_FILE})

source_group ("Generated Files" FILES ${MOC_SOURCES} ${QRC_SOURCES} ${QM_FILE})
source_group ("Translation Files" FILES ${TS_FILE})
add_executable(${PROJECT_NAME}
        ${SOURCE_FILE_LIST}
        ${HEADER_FILES}
        ${MOC_SOURCES}
        ${QRC_SOURCES}
        ${TS_FILE}
        ${QM_FILE}
        ${WINRC_FILES}
        )
if (UNIX)
    set_target_properties (${PROJECT_NAME} PROPERTIES
        COMPILE_FLAGS -fPIE
        )
endif(UNIX)

if(WIN32)
  set_property(TARGET ${PROJECT_NAME} PROPERTY WIN32_EXECUTABLE true)
endif(WIN32)
