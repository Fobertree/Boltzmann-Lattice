#version 330 core

// colormap glsl shader
in vec2 TexCoords; // for sampler2D
out vec4 FragColor;

uniform sampler2D dataTexture; // 2D matrix data as texture
uniform float minValue;
uniform float maxValue;

vec3 applyColormap(float value) {
    float normalizedValue = (value - minValue) / (maxValue - minValue);

    normalizedValue = clamp(normalizedValue, 0.0, 1.0);

    // default color pallete
    vec3 blue = vec3(0.0, 0.0, 1.0);
    vec3 white = vec3(1.0, 1.0, 1.0);
    vec3 red = vec3(1.0, 0.0, 0.0);

    // interpolation
    if (normalizedValue < 0.5) {
        // Interpolate between blue and white
        return mix(blue, white, normalizedValue * 2.0);
    } else {
        // Interpolate between white and red
        return mix(white, red, (normalizedValue - 0.5) * 2.0);
    }
}

void main() {
    float value = texture(dataTexture, TexCoords).r; // sample data
    FragColor = vec4(applyColormap(value), 1.0);
}