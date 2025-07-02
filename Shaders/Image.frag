uniform sampler2D iChannel0;
uniform sampler2D iChannel1;
uniform sampler2D iChannel2;
uniform sampler2D iChannel3;
uniform float iTime;
uniform vec4 iMouse;
uniform vec2 iResolution;
out vec4 C;


mat3 edge = mat3(-1.,-1.,-1.,
                                     -1.,8.,-1.,
                                     -1.,-1.,-1.);

mat3 gauss = mat3(1.,1.,1.,
                                     1.,1.,1.,
                                     1.,1.,1.);

mat3 laplacian = mat3(0.,-1.,0.,
                                     -1.,4.,-1.,
                                     0.,-1.,0.);

mat3 sobelH = mat3(1.,0.,-1.,
                                     2.,0.,-2.,
                                     1.,0.,-1.);

mat3 sobelV = mat3(1.,2.,1.,
                                     0.,0.,0.,
                                     -1.,-2.,-1.);



float Filter33(vec2 pos,mat3 kernel)
{
         float r=0.0;
                 for(int i=-1;i<2;i++)
                    for (int j =-1;j<2;j++)
                    {
                        vec3 c =texture(iChannel0,(pos+vec2(i,j))/iResolution.xy).xyz;
                        float tmp = dot(c,vec3(0.177,0.812,0.0106));
                        r+=tmp*kernel[i+1][j+1];
                    }
         return r/9.;
}


vec3 Filter33c(vec2 pos,mat3 kernel)
{
         vec3 r=vec3(0.0);
                 for(int i=-1;i<2;i++)
                    for (int j =-1;j<2;j++)
                    {
                        vec3 c =texture(iChannel0,(pos+vec2(i,j))/iResolution.xy).xyz;

                        r+=c*kernel[i+1][j+1];
                    }
         return r/9.;
}



void main() {
    vec2 uv = gl_FragCoord.xy / iResolution.xy;

    vec3 baseColor = texture(iChannel0, uv).rgb;

    // --- Bloom inline refinado ---
    vec3 bloom = vec3(0.0);
    float threshold = 3.0;       // só ativa em regiões realmente brilhantes
    float strength = 0.3;        // força discreta
    float radius = 1.5;
    float totalWeight = 0.0;

    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            vec2 offset = vec2(float(x), float(y)) / iResolution * radius;
            vec3 s = texture(iChannel0, uv + offset).rgb;

            float brightness = max(max(s.r, s.g), s.b);
            if (brightness > threshold) {
                float dist = length(vec2(x, y));
                float weight = exp(-dist * dist * 1.2); // peso gaussiano leve
                bloom += s * weight;
                totalWeight += weight;
            }
        }
    }

    if (totalWeight > 0.0) {
        bloom /= totalWeight;
        baseColor += bloom * strength;
    }

    // --- Pós-processamento ---
    // baseColor = baseColor / (baseColor + vec3(1.0)); // tone mapping
    // baseColor = pow(baseColor, vec3(1.0 / 2.2));      // gamma correction

    C = vec4(baseColor, 1.0);
}

