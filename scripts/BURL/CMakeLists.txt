cmake_minimum_required(VERSION 3.20)
project(BURL CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

# ── Collect sources ───────────────────────────────────────────────────────────
file(GLOB SOURCES "BURL/*.cpp")

# ── Define the executable ─────────────────────────────────────────────────────
add_executable(burl ${SOURCES})

# ── Include directories ───────────────────────────────────────────────────────
target_include_directories(burl PRIVATE
    ${CMAKE_SOURCE_DIR}/BURL
    ${CMAKE_SOURCE_DIR}/BURL/Eigen
)

# ── Compiler flags ────────────────────────────────────────────────────────────
target_compile_options(burl PRIVATE
    -O3
    -Wall
    -Wextra
)