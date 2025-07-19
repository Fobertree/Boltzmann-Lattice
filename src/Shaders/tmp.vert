#version 330 core
// skip VBO for full-screen shader

layout(location = 0) in vec2 vertexIn;
out vec2 textureCoord;

void main() {
    // [-1,1] -> [0, -1]
    textureCoord = vertexIn * 0.5 + 0.5;
    gl_Position = vec4(vertexIn, 0.0, 1.0);
}