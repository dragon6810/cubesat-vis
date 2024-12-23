# Compiling

## Debian
    sudo apt update
    sudo apt install libglfw3-dev libglew-dev
    make

## MacOS (homebrew)
    brew install glfw glew
    sudo ln -s /opt/homebrew/Cellar/glew/2.2.0_1/include/GL /usr/local/include/GL
    sudo ln -s /opt/homebrew/Cellar/glfw/3.4/include/GLFW /usr/local/include/GLFW
    sudo ln -sf /opt/homebrew/Cellar/glew/2.2.0_1/lib/libGLEW.dylib /usr/local/lib/libGLEW.dylib
    sudo ln -sf /opt/homebrew/Cellar/glfw/3.4/lib/libglfw.dylib /usr/local/lib/libglfw.dylib
    sudo ln -sf /opt/homebrew/Cellar/glfw/3.4/lib/libglfw3.a /usr/local/lib/libglfw3.a
    make

replase 2.2.0_1 with whatever version you have, check with
    brew info glew
same goes for glfw