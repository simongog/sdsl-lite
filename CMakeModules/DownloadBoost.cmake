ExternalProject_Add(
    boost_external
    GIT_REPOSITORY https://github.com/boostorg/boost.git
    GIT_TAG boost-1.59.0
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ./bootstrap.sh
    BUILD_COMMAND ./b2 --prefix=<INSTALL_DIR> --without-python install
    #BUILD_COMMAND ./b2 --show-libraries
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(boost_external source_dir install_dir)

file(MAKE_DIRECTORY "${install_dir}/include")

set(${package_found_prefix}_CMAKE_DEP boost_external)
set(${package_found_prefix}_LIBRARIES
    "${install_dir}/lib/libboost_system.a"
    "${install_dir}/lib/libboost_filesystem.a"
    "${install_dir}/lib/libboost_program_options.a"
)
set(${package_found_prefix}_INCLUDE_DIRS "${install_dir}/include")
