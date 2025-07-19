#version 330

#ifdef GL_ES
precision mediump float;
#endif

out vec4 color;

void main() {
    vec3 rgb = (1., 0., 0.);
    color = vec4(rgb,1.);
}