cmake_minimum_required(VERSION 3.20)
project(BURLc CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# ── Collect sources ───────────────────────────────────────────────────────────
file(GLOB SOURCES "BURL_coverage/*.cpp")

# ── Define the executable ─────────────────────────────────────────────────────
add_executable(burlc ${SOURCES})

# ── Include directories ───────────────────────────────────────────────────────
target_include_directories(burlc PRIVATE
    ${CMAKE_SOURCE_DIR}/BURL_coverage
    ${CMAKE_SOURCE_DIR}/BURL_coverage/Eigen
)

# ── Compiler flags ────────────────────────────────────────────────────────────
target_compile_options(burlc PRIVATE
    -O3
    -Wall
    -Wextra
)
