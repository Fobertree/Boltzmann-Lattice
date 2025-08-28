#version 330 core

// colormap glsl shader
in vec2 TexCoords; // for sampler2D
out vec4 FragColor;

uniform sampler2D dataTexture; // 2D matrix data as texture
uniform float minValue;
uniform float maxValue;

vec3 applyColormap(float value) {
    float normalizedValue = (value - minValue) / (maxValue - minValue);
    vec3 hexColor1 = vec3(0.69, 0.576, 0.929); // 176, 147, 237
    vec3 hexColor2 = vec3(0.,0.,0.);
    // interpolation
    // this shader can be done better
    return mix(hexColor2, hexColor1, value * 2.3);
}

void main() {
    float value = texture(dataTexture, TexCoords).r; // sample data
    FragColor = vec4(applyColormap(value), 1.0);
}