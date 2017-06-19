cmake_minimum_required (VERSION 2.8.3) 
project (vibro-transport) 

include (../qt5_executable.cmake)

add_dependencies (${PROJECT_NAME} ode_num_int)

target_link_libraries (${PROJECT_NAME}
    ode_num_int
    ${Qt5Widgets_LIBRARIES}
    ${Qt5Core_LIBRARIES}
    )

set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY
  COMPILE_DEFINITIONS $<$<CONFIG:Debug>:_DEBUG>
)
