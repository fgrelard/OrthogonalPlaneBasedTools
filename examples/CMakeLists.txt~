SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options system filesystem)
if ( Boost_FOUND )
    message(STATUS "Boost and boost_program_options found.")
    include_directories(${Boost_INCLUDE_DIRS})
    SET(MyDGtalContrib ${DGtalToolsLibDependencies}
            ${Boost_LIBRAIRIES}
            ${Boost_PROGRAM_OPTIONS_LIBRARY}
            ${Boost_FILESYSTEM_LIBRARY}
            ${Boost_SYSTEM_LIBRARY})
    SET(MyDGtalContrib ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )

add_subdirectory(vcm)
add_subdirectory(converters)
add_subdirectory(geometry)
add_subdirectory(hessian)