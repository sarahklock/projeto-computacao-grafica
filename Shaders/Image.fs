#version 330 core

uniform sampler2D iChannel0;
uniform vec2 iResolution;
out vec4 C;

void main() {
    vec2 uv = gl_FragCoord.xy / iResolution.xy;

    // Lê a cor do framebuffer anterior (BufferA, por exemplo)
    vec3 color = texture(iChannel0, uv).rgb;

    // Tone mapping de Reinhard
    color = color / (color + vec3(1.0));

    // Correção gama (opcional, melhora aparência)
    color = pow(color, vec3(1.0 / 2.2));

    C = vec4(color, 1.0);
}
