SET(sdsl_tests "")
SET(generated_files "")

FILE(STRINGS ${CMAKE_CURRENT_LIST_DIR}/download.config download_config_lines REGEX "^[^#].*")
FOREACH(download_config_line ${download_config_lines})
    LIST(GET download_config_line 0 download_file)
    LIST(GET download_config_line 1 download_url)
    STRING(REPLACE ".tar.gz" "" download_file_we ${download_file})
    SET(abs_test_case_name ${CMAKE_CURRENT_LIST_DIR}/test_cases/${download_file_we})
    
    ADD_CUSTOM_COMMAND(OUTPUT ${abs_test_case_name}
                      COMMAND ${CMAKE_COMMAND} -Ddownload_url=${download_url} -Dtest_case_dir=${CMAKE_CURRENT_LIST_DIR}/test_cases -P ${CMAKE_CURRENT_LIST_DIR}/download.cmake
                      WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
                      COMMENT "Checking or downloading test case ${download_file}."
                      VERBATIM
    )
    ADD_CUSTOM_TARGET(${download_file_we}
                      DEPENDS ${abs_test_case_name})
    LIST(APPEND generated_files ${abs_test_case_name})
ENDFOREACH(download_config_line)


FILE(GLOB SDSL_HEADER_ALL ${CMAKE_CURRENT_BINARY_DIR}/../include/sdsl/*.hpp)
STRING(REGEX REPLACE "([^;]+)[;]" "#include \"\\1\"\\n" SDSL_INCLUDE_ALL "${SDSL_HEADER_ALL};")                                                                       
STRING(REGEX REPLACE "${CMAKE_CURRENT_BINARY_DIR}/\\.\\./include/" "" SDSL_INCLUDE_ALL "${SDSL_INCLUDE_ALL}")                                                                       

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/compile_test.cpp.cmake"
               "${CMAKE_CURRENT_SOURCE_DIR}/compile_test.cpp" @ONLY)

FILE(GLOB test_sources ${CMAKE_CURRENT_SOURCE_DIR}/*_test.cpp)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../include
                    ${CMAKE_CURRENT_BINARY_DIR}/../external/libdivsufsort/include
                    ${gtest_SOURCE_DIR}/include
)

LINK_LIBRARIES(sdsl gtest divsufsort divsufsort64)

ADD_EXECUTABLE(int_vector_generator EXCLUDE_FROM_ALL int_vector_generator.cpp)
ADD_EXECUTABLE(bit_vector_generator EXCLUDE_FROM_ALL bit_vector_generator.cpp)
ADD_EXECUTABLE(replace_int_vector_value EXCLUDE_FROM_ALL replace_int_vector_value.cpp)

SET(tmp_dir ${CMAKE_CURRENT_LIST_DIR}/tmp)
FOREACH(test_source ${test_sources})
    GET_FILENAME_COMPONENT(test_name_we ${test_source} NAME_WE)

    STRING(REGEX REPLACE "_" "-" test_name_hy ${test_name_we})

    ADD_EXECUTABLE(${test_name_we} EXCLUDE_FROM_ALL ${test_name_we}.cpp)

    SET(config_file ${CMAKE_CURRENT_SOURCE_DIR}/${test_name_we}.config)
# (1) Handle tests without .config file
    IF(NOT EXISTS ${config_file})
        ADD_CUSTOM_TARGET(${test_name_hy}
            COMMAND $<TARGET_FILE:${test_name_we}> ${tmp_dir}
            DEPENDS test/${test_name_we} ${test_name_we}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Execute ${test_name_hy}."
        )
        LIST(APPEND sdsl_tests ${test_name_hy})
# (2) Handle tests with .config file        
    ELSE(NOT EXISTS ${config_file})
        SET(test_targets "")
        FILE(STRINGS ${config_file} config_lines REGEX "^[^#].+")
        SET(config_line_cnt 0)
        FOREACH(config_line ${config_lines})
            MATH(EXPR config_line_cnt ${config_line_cnt}+1)
# (2.1) Handle tests with simple .config file (only one input file)
            IF("${config_line}" MATCHES "^[^;]+$")
                GET_FILENAME_COMPONENT(test_case_name ${config_line} NAME)
                SET(abs_test_case_name ${CMAKE_CURRENT_LIST_DIR}/test_cases/${test_case_name})
               
                SET(int_vec_regex "int-vec\\.[0-9]+\\.[0-9]+\\..+")
                SET(bit_vec_regex "bit-vec\\.*")
                IF("${test_case_name}" MATCHES "${int_vec_regex}")
                    STRING(REPLACE "int-vec." "" test_case_suf ${test_case_name})
                    #  insert semicolons to get a list :)
                    STRING(REPLACE "." ";" test_case_param ${test_case_suf}) 

                    ADD_CUSTOM_COMMAND(OUTPUT ${abs_test_case_name}
                                       COMMAND $<TARGET_FILE:int_vector_generator> ${abs_test_case_name}  ${test_case_param} 
                                       DEPENDS int_vector_generator
                                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                                       COMMENT "Generating test case ${test_case_name}."
                                       VERBATIM
                    )
                    LIST(APPEND generated_files ${abs_test_case_name})
                ELSEIF("${test_case_name}" MATCHES "${bit_vec_regex}")   
                    STRING(REPLACE "bit-vec." "" test_case_param ${test_case_name})
                    ADD_CUSTOM_COMMAND(OUTPUT ${abs_test_case_name}
                                       COMMAND $<TARGET_FILE:bit_vector_generator> ${abs_test_case_name} ${test_case_param} 
                                       DEPENDS bit_vector_generator
                                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                                       COMMENT "Generating test case ${test_case_name}."
                                       VERBATIM
                    )
                    LIST(APPEND generated_files ${abs_test_case_name})
                ENDIF("${test_case_name}" MATCHES "${int_vec_regex}")

                SET(test_target ${test_name_hy}_${test_case_name})
                ADD_CUSTOM_TARGET(${test_target}
                    COMMAND $<TARGET_FILE:${test_name_we}> ${abs_test_case_name} ${tmp_dir}/${test_name_we}_${test_case_name} ${tmp_dir}
                    DEPENDS ${abs_test_case_name} ${tmp_dir} ${test_name_we}
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR} 
                    COMMENT "Execute ${test_name_hy} on ${test_case_name}."
                )
                LIST(APPEND test_targets ${test_target})
            ELSE("${config_line}" MATCHES "^[^;]+$")
#               Handle Cs[a|t]IntTest
                IF ("${test_name_we}" MATCHES "cs[a|t]_int_test") 
                    LIST(GET config_line 0 test_case_name)
                    SET(abs_test_case_name ${CMAKE_CURRENT_LIST_DIR}/test_cases/${test_case_name})

                    SET(int_vec_sa_regex "int-vec-sa\\.[0-9]+\\.[0-9]+\\..+")
                    IF("${test_case_name}" MATCHES "${int_vec_sa_regex}")
                        STRING(REPLACE "int-vec-sa." "" test_case_suf ${test_case_name})
                        #  insert semicolons to get a list :)
                        STRING(REPLACE "." ";" test_case_param ${test_case_suf}) 

                        ADD_CUSTOM_COMMAND(OUTPUT ${abs_test_case_name}
                                           COMMAND $<TARGET_FILE:int_vector_generator> ${abs_test_case_name}  ${test_case_param} 
                                           COMMAND $<TARGET_FILE:replace_int_vector_value> ${abs_test_case_name} 0 1
                                           DEPENDS int_vector_generator replace_int_vector_value 
                                           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                                           COMMENT "Generating test case ${test_case_name}."
                                           VERBATIM
                        )
                        LIST(APPEND generated_files ${abs_test_case_name})
                    ENDIF("${test_case_name}" MATCHES "${int_vec_sa_regex}")


                    LIST(GET config_line 1 num_byte)
                    SET(test_target ${test_name_hy}_${test_case_name}_${num_byte})
                    ADD_CUSTOM_TARGET(${test_target}
                        COMMAND $<TARGET_FILE:${test_name_we}> ${abs_test_case_name} ${num_byte} ${tmp_dir}/${test_name_we}_${tase_case_name}_${num_byte} ${tmp_dir}
                        DEPENDS ${abs_test_case_name} ${tmp_dir} ${test_name_we}
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                        COMMENT "Execute ${test_name_hy} on ${test_case_name}."        
                    )
                    LIST(APPEND test_targets ${test_target})
#               Handle K2TreapTest    
                ELSEIF("${test_name_we}" MATCHES "k2_treap_test")
                    SET(concat_config "")
#                    STRING(CONCAT concat_config ${config_line})
                    FOREACH(item ${config_line})
                        SET(concat_config "${concat_config}${item}")
                    ENDFOREACH(item)
                    SET(k2dim x y w)
                    LIST(GET config_line 0 test_case_name)
                    SET(test_case_name "k2-${concat_config}")
                    SET(abs_test_case_name ${CMAKE_CURRENT_LIST_DIR}/test_cases/${test_case_name})
                    SET(abs_test_case_name_list "")
                    FOREACH(i 0 1 2)
                        LIST(GET k2dim ${i} dim)
                        SET(local_abs_test_case_name ${abs_test_case_name}.${dim})
                        LIST(APPEND abs_test_case_name_list ${local_abs_test_case_name})
                        LIST(GET config_line ${i} test_case_suf)
                        STRING(REPLACE "." ";" test_case_param ${test_case_suf})
                        ADD_CUSTOM_COMMAND(OUTPUT ${local_abs_test_case_name}
                                COMMAND $<TARGET_FILE:int_vector_generator> ${local_abs_test_case_name} ${test_case_param}
                                DEPENDS int_vector_generator
                                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                                COMMENT "Generating test case ${test_case_name}."
                                VERBATIM        
                        )
                    ENDFOREACH(i)
                    ADD_CUSTOM_TARGET(${test_case_name} 
                            DEPENDS ${abs_test_case_name_list})
                    SET(test_target ${test_name_hy}_${test_case_name})
                    ADD_CUSTOM_TARGET(${test_target}
                                COMMAND $<TARGET_FILE:${test_name_we}> ${abs_test_case_name} ${tmp_dir}/${test_name_we}_${test_case_name}
                                DEPENDS ${test_case_name} ${tmp_dir} ${test_name_we}
                                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                                COMMENT "Execute ${test_name_hy} on ${test_case_name}."
                            )
                    LIST(APPEND test_targets ${test_target})
                ENDIF("${test_name_we}" MATCHES "cs[a|t]_int_test")
            ENDIF("${config_line}" MATCHES "^[^;]+$")
        ENDFOREACH(config_line)
        ADD_CUSTOM_TARGET(${test_name_hy} DEPENDS ${test_targets}) 
        LIST(APPEND sdsl_tests ${test_name_hy})
    ENDIF(NOT EXISTS ${config_file})
ENDFOREACH(test_source)

ADD_CUSTOM_TARGET(test-sdsl DEPENDS ${sdsl_tests})

ADD_CUSTOM_TARGET(clean-test
                 COMMAND ${CMAKE_COMMAND} -E remove -f ${generated_files}
                 COMMENT "Remove generated and downloaded test inputs.")


# TODO: special handling for in-memory test
