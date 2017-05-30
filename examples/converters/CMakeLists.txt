if (WITH_ITK)
  SET(DGTAL_CONTRIB_EXAMPLES
    dicom2vol
    vol2dicom
    )
endif(WITH_ITK)

FOREACH(FILE ${DGTAL_CONTRIB_EXAMPLES})
  add_executable(${FILE} ${FILE})
  target_link_libraries (${FILE} DGtal DGtalIO ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
  set_target_properties( ${FILE} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${EXAMPLE_OUTPUT_PATH} )
ENDFOREACH(FILE)