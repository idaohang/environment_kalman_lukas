set ( SOURCES
        "src/TestClass.cpp"
        )

set (HEADERS
        "include/filter/system_model.h"
        "include/filter/measurement_model.h"
        "include/TestClass.h")

include_directories(
        "include"

        # Eigen3 include paths
        "/usr/include/eigen3"
        "/usr/local/include/eigen3"

        # Kalman Lib paths
        "/usr/include"
        "/usr/local/include"
)

add_library (environment_kalman_lukas MODULE ${SOURCES} ${HEADERS})
target_link_libraries(environment_kalman_lukas PRIVATE lmscore math_lib sensor_utils)
