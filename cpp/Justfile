default: debug

build_dir := "build"

debug:
    @if [ -d "{{build_dir}}" ]; then rm -rf {{build_dir}}; fi
    @cmake -B {{build_dir}} -G Ninja
    @cmake --build {{build_dir}}

release:
    @if [ -d "{{build_dir}}" ]; then rm -rf {{build_dir}}; fi
    @cmake -B {{build_dir}} -G Ninja -DCMAKE_BUILD_TYPE=Release
    @cmake --build {{build_dir}}
