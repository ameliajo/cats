
add_executable( gasket.helperFiles.test
                acon2.test.cpp
                bessl.test.cpp
                ftrans.test.cpp
                gtg.test.cpp
                rconv.test.cpp
                scint.test.cpp
                sterp.test.cpp )
target_compile_options( gasket.helperFiles.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${cats_appended_flags} )
target_link_libraries( gasket.helperFiles.test PUBLIC cats ) 
add_test( NAME gasket.helperFiles COMMAND gasket.helperFiles.test )