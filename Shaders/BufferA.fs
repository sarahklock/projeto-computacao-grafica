in vec3 aColor;
in vec4 aPosition;
out vec4 C;
uniform sampler2D iChannel0;
uniform sampler2D iChannel1;
uniform sampler2D iChannel2;
uniform sampler2D iChannel3;
uniform vec2 iResolution;
uniform vec4 iMouse;
uniform float iTime;
uniform int iFrame;



#define MAX_STEPS 100
#define MAX_DIST 200.
#define SURF_DIST .01
#define EPSILON .01
#define PI 3.14159265359
float dot2( in vec2 v ) { return dot(v,v); }
float dot2( in vec3 v ) { return dot(v,v); }
float ndot( in vec2 a, in vec2 b ) { return a.x*b.x - a.y*b.y; }

struct Surface {
    float distanceToSurface; // Distância do raio até o ponto de interseção 
    vec3 surfaceColor;       // Cor da superfície no ponto de interseção
    float ambientCoeff;      // Coeficiente de iluminação ambiente (Ka)
    float diffuseCoeff;      // Coeficiente de iluminação difusa (Kd)
    float specularCoeff;     // Coeficiente de iluminação especular (Ks)
    int objectId;            // ID do objeto atingido 
};


vec3 getSkyColor(vec3 rayDirection, vec3 sunPosition, float sunAngle) {
    
    // Normaliza a altura do sol para [0,1]: -3 → 0, 0 → 0.5, +3 → 1
    float sunHeight = clamp(sunPosition.y / 3.0 * 0.5 + 0.5, 0.0, 1.0);

    vec3 nightColor = vec3(0.02, 0.02, 0.08);  // noite
    vec3 dawnColor  = vec3(0.8, 0.4, 0.2);     // amanhecer/entardecer
    vec3 dayColor   = vec3(0.4, 0.7, 1.0);     // dia

    // Se sunHeight está entre 0 e 0.5, interpola entre noite e amanhecer
    // Se entre 0.5 e 1.0, interpola entre amanhecer e dia

    vec3 skyColor;
    if (sunHeight < 0.5) {
        float t = smoothstep(0.0, 0.5, sunHeight);
        skyColor = mix(nightColor, dawnColor, t);
    } else {
        float t = smoothstep(0.5, 1.0, sunHeight);
        skyColor = mix(dawnColor, dayColor, t);
    }

    // Atenuação perto do horizonte para escurecer
    float verticalFactor = 1.0 / max(rayDirection.y, 0.01);
    vec3 attenuation = exp2(-verticalFactor * vec3(0.4, 0.6, 1.0));
    skyColor = mix(skyColor, vec3(0.0), attenuation);

    return skyColor;
}


const vec3 COLOR_BACKGROUND = vec3(0.25, 0.1, 0.15);

mat4 translateMatrix(vec3 translation) {
    mat4 result = mat4(
        vec4(1.0, 0.0, 0.0, 0.0),
        vec4(0.0, 1.0, 0.0, 0.0),
        vec4(0.0, 0.0, 1.0, 0.0),
        vec4(translation.x, translation.y, translation.z, 1.0)
    );
    return result;
}


// vertical
float signedDistanceCylinderVertical(vec3 point, vec2 halfSize) {
    // Calcula a distância em xz e y para o cilindro com meio tamanho halfSize
    vec2 distanceVec = abs(vec2(length(point.xz), point.y)) - halfSize;
    return min(max(distanceVec.x, distanceVec.y), 0.0) + length(max(distanceVec, 0.0));
}


Surface unionSurfaces(Surface surface1, Surface surface2) {
    if (surface1.distanceToSurface <= surface2.distanceToSurface)
        return surface1;
    else
        return surface2;
}


Surface getDist(vec3 point) {
    // Define o cilindro vertical
    vec2 cylinderSize = vec2(0.35, 1.2); // raio = 0.35, altura = 1.2
    vec3 cylinderCenter = vec3(0.0, 1.2, 4.0); // posição central do cilindro

    Surface cylinderSurface;
    cylinderSurface.distanceToSurface = signedDistanceCylinderVertical(point - cylinderCenter, cylinderSize);
    cylinderSurface.surfaceColor = vec3(1.0); // cor padrão branca, será sobrescrita pela textura
    cylinderSurface.ambientCoeff = 0.2;
    cylinderSurface.diffuseCoeff = 0.4;
    cylinderSurface.specularCoeff = 0.4;
    cylinderSurface.objectId = 10;  // ID exclusivo para o cilindro (texturizado)

    // Define o chão
    Surface groundSurface;
    groundSurface.distanceToSurface = point.y; // plano horizontal no y=0
    groundSurface.surfaceColor = vec3(0.2, 0.8, 0.2); // verde
    groundSurface.ambientCoeff = 0.1;
    groundSurface.diffuseCoeff = 0.5;
    groundSurface.specularCoeff = 0.1;
    groundSurface.objectId = 20; // ID exclusivo para o chão

    // Combina cilindro e chão usando união booleana
    Surface combinedSurface = unionSurfaces(groundSurface, cylinderSurface);

    // Parâmetros para o sol (esfera animada)
    float sunRadius = 0.3;
    float orbitRadius = 6.0; // raio da orbita do sol e da lua
    float sunOrbitAngle = iTime * 0.2; // velocidade angular do sol
    float orbitCenter = 1.0; //define a altura do sol e da lua

    vec3 sunPosition = vec3(
        orbitRadius * cos(sunOrbitAngle), 
        3.0 * sin(sunOrbitAngle) + orbitCenter, 
        4.0
    );

    // Define a superfície do sol (esfera)
    Surface sunSurface;
    sunSurface.distanceToSurface = length(point - sunPosition) - sunRadius;
    sunSurface.surfaceColor = vec3(1.0, 0.9, 0.5); // tom amarelado
    sunSurface.ambientCoeff = 0.0;
    sunSurface.diffuseCoeff = 0.0;
    sunSurface.specularCoeff = 0.0;
    sunSurface.objectId = 30; // ID exclusivo para o sol

    // Combina o sol com a cena anterior
    combinedSurface = unionSurfaces(sunSurface, combinedSurface);

    // Parâmetros para a lua
    float moonRadius = 0.3;
    vec3 moonPosition = vec3(
        -orbitRadius * cos(sunOrbitAngle),
        -3.0 * sin(sunOrbitAngle) + orbitCenter,
        4.0
    );

    // Define a superfície da lua (esfera oposta ao Sol)
    Surface moonSurface;
    moonSurface.distanceToSurface = length(point - moonPosition) - moonRadius;
    moonSurface.surfaceColor = vec3(0.6, 0.7, 1.0); // tom azulado
    moonSurface.ambientCoeff = 0.1;
    moonSurface.diffuseCoeff = 0.5;
    moonSurface.specularCoeff = 0.2;
    moonSurface.objectId = 40;

    // Combina lua com a cena
    combinedSurface = unionSurfaces(moonSurface, combinedSurface);

    return combinedSurface;
}

Surface rayMarching(vec3 cameraPosition, vec3 rayDirection) {
    float totalDistance = 0.0;
    vec3 currentPosition = cameraPosition + totalDistance * rayDirection;

    // Avalia a distância até o objeto mais próximo na cena
    Surface surfaceHit = getDist(currentPosition);

    // Posição animada do sol (deve bater com a usada em getDist)
    float sunAngle = iTime * 0.2;
    float orbitRadius = 6.0;
    float lightY = 3.0 * sin(sunAngle);

    vec3 sunPosition  = vec3( orbitRadius * cos(sunAngle),  lightY, 4.0);
    vec3 moonPosition = vec3(-orbitRadius * cos(sunAngle), -lightY, 4.0);

    float sunGlow = 0.0;
    float moonGlow = 0.0;

    int stepCount = 0;

    // Loop de ray marching: avança até encontrar a superfície ou exceder os limites
    while ((stepCount < MAX_STEPS) && (surfaceHit.distanceToSurface > SURF_DIST) && (totalDistance < MAX_DIST)) {
        float distanceToSun = length(currentPosition - sunPosition);
        sunGlow += exp(-distanceToSun * 4.0) * 0.08;


        float distanceToMoon = length(currentPosition - moonPosition);
        moonGlow += exp(-distanceToMoon * 10.0) * 0.01;

        totalDistance += surfaceHit.distanceToSurface;
        currentPosition = cameraPosition + totalDistance * rayDirection;
        surfaceHit = getDist(currentPosition);
        stepCount++;
    }

    // Se o raio não colidiu com nenhum objeto (céu ou fundo)
    if ((stepCount > MAX_STEPS) || (totalDistance > MAX_DIST)) {
        surfaceHit.surfaceColor = getSkyColor(rayDirection, sunPosition, sunAngle)
    + sunGlow * vec3(1.0, 0.9, 0.6) * 1.5
    + moonGlow * vec3(0.6, 0.7, 1.0);

        surfaceHit.distanceToSurface = MAX_DIST;
    } else {
        surfaceHit.distanceToSurface = totalDistance;

        // Caso o raio tenha atingido diretamente o sol (ID 30)
        if (surfaceHit.objectId == 30) {
            surfaceHit.surfaceColor = vec3(1.0, 0.9, 0.6) * 8.0; // brilho intenso do sol visível
        }
    }

    return surfaceHit;
}

vec3 estimateNormal(vec3 position) {
    float baseDistance = getDist(position).distanceToSurface;

    float deltaX = getDist(vec3(position.x + EPSILON, position.y, position.z)).distanceToSurface - baseDistance;
    float deltaY = getDist(vec3(position.x, position.y + EPSILON, position.z)).distanceToSurface - baseDistance;
    float deltaZ = getDist(vec3(position.x - EPSILON, position.y, position.z + EPSILON)).distanceToSurface - baseDistance;

    vec3 gradient = vec3(deltaX, deltaY, deltaZ);

    return normalize(gradient);
}

mat3 setCamera(vec3 cameraPosition, vec3 targetPosition)
{
    // Vetor direção da câmera (do ponto da câmera para o alvo)
    vec3 cameraDirection = normalize(targetPosition - cameraPosition);

    // Vetor "direita" da câmera, perpendicular ao vetor direção e ao vetor up global (0,1,0)
    vec3 cameraRight = cross(cameraDirection, vec3(0.0, 1.0, 0.0));

    // Vetor "up" da câmera, perpendicular aos vetores direita e direção
    vec3 cameraUp = cross(cameraRight, cameraDirection);

    // Monta a matriz de orientação da câmera:
    // Note que o vetor direita é invertido para adequar o sistema de coordenadas
    return mat3(-cameraRight, cameraUp, cameraDirection);
}


float calcSoftshadow(
    in vec3 rayOrigin, 
    in vec3 rayDirection, 
    in float minDistance, 
    in float maxDistance, 
    int ignoredObjectId
) {
    // Cálculo para limitar o alcance da sombra com base no plano do chão (bounding volume)
    float floorIntersection = (maxDistance * 0.75 - rayOrigin.y) / rayDirection.y;
    if (floorIntersection > 0.0)
        maxDistance = min(maxDistance, floorIntersection);

    float shadowFactor = 1.0;
    float currentDistance = minDistance;

    for (int step = 0; step < 24; step++) {
        vec3 samplePoint = rayOrigin + currentDistance * rayDirection;
        Surface surfaceSample = getDist(samplePoint);

        // Ignora objetos com ID específico (ex: o próprio Sol)
        if (surfaceSample.objectId == ignoredObjectId) {
            currentDistance += 0.05;
            continue;
        }

        float distanceToSurface = surfaceSample.distanceToSurface;
        float attenuation = clamp(8.0 * distanceToSurface / currentDistance, 0.0, 1.0);
        shadowFactor = min(shadowFactor, attenuation);

        currentDistance += clamp(distanceToSurface, 0.01, 0.25);

        // Se sombra for quase total ou ultrapassamos a distância máxima, encerramos
        if (shadowFactor < 0.001 || currentDistance > maxDistance)
            break;
    }

    // Suavização final com interpolação cúbica (Hermite)
    shadowFactor = clamp(shadowFactor, 0.0, 1.0);
    return shadowFactor * shadowFactor * (3.0 - 2.0 * shadowFactor);
}

vec3 getLight(vec3 surfacePosition, Surface surfaceData, vec3 cameraPosition) {
    // Se não há superfície (raymarching não encontrou nada), retorna a cor do fundo
    if (surfaceData.distanceToSurface == MAX_DIST)
        return surfaceData.surfaceColor;

    // Posição e cor do Sol e da Lua (animado)
    float sunAngle = iTime * 0.2;
    float orbitRadius = 6.0;
    float lightY = 3.0 * sin(sunAngle);
    vec3 sunPosition  = vec3( orbitRadius * cos(sunAngle),  lightY, 4.0);
    vec3 moonPosition = vec3(-orbitRadius * cos(sunAngle), -lightY, 4.0);

    vec3 sunColor = vec3(1.0, 0.9, 0.6);
    vec3 moonColor = vec3(0.5, 0.6, 1.0); // luz azulada fraca

    // Segunda fonte de luz (vermelha, fixa)
    vec3 redLightColor = vec3(1.0, 0.0, 0.0);
    vec3 redLightPosition = vec3(-3.0, 4.0, 4.0);

    // Direções das luzes para o ponto
    vec3 sunDirection = normalize(sunPosition - surfacePosition);
    vec3 redLightDirection = normalize(redLightPosition - surfacePosition);
    vec3 moonDirection = normalize(moonPosition - surfacePosition);
    vec3 viewDirection = normalize(cameraPosition - surfacePosition);

    // Normal estimada no ponto da superfície
    vec3 surfaceNormal = estimateNormal(surfacePosition);

    // Vetores de reflexão especular
    vec3 sunReflection = normalize(reflect(-sunDirection, surfaceNormal));
    vec3 redReflection = normalize(reflect(-redLightDirection, surfaceNormal));
    vec3 moonReflection = normalize(reflect(-moonDirection, surfaceNormal));

    // Iluminação difusa (dot produto entre luz e normal)
    float sunDiffuse = clamp(dot(surfaceNormal, sunDirection), 0.0, 1.0);
    float redDiffuse = clamp(dot(surfaceNormal, redLightDirection), 0.0, 1.0);
    float moonDiffuse = clamp(dot(surfaceNormal, moonDirection), 0.0, 1.0);

    // Sombras suaves (evita sombra sobre o próprio sol)
    float sunShadowFactor = 1.0;
    float redShadowFactor = 1.0;
    // if (surfaceData.objectId != 30) {
    //     sunShadowFactor = calcSoftshadow(surfacePosition + 10.0 * EPSILON * sunDirection, sunDirection, 0.1, 3.0, 30);
    //     redShadowFactor = calcSoftshadow(surfacePosition + 10.0 * EPSILON * redLightDirection, redLightDirection, 0.1, 3.0, -1);
    // }

    // Reflexão especular (Phong)
    float shininess = 7.0;
    vec3 sunSpecular = vec3(0.0);
    vec3 redSpecular = vec3(0.0);

    float sunSpecDot = dot(sunReflection, viewDirection);
    if (sunSpecDot > 0.0)
        sunSpecular = sunColor * surfaceData.specularCoeff * pow(sunSpecDot, shininess);

    float redSpecDot = dot(redReflection, viewDirection);
    if (redSpecDot > 0.0)
        redSpecular = redLightColor * surfaceData.specularCoeff * pow(redSpecDot, shininess);

    float moonSpecDot = dot(moonReflection, viewDirection);
    vec3 moonSpecular = vec3(0.0);
    if (moonSpecDot > 0.0)
        moonSpecular = moonColor * surfaceData.specularCoeff * pow(moonSpecDot, shininess);

    // Textura cilíndrica (aplicada ao objeto com ID 10)
    if (surfaceData.objectId == 10) {
        vec3 localPosition = surfacePosition - vec3(0.0, 1.2, 4.0); // usa o centro real do cilindro
        float theta = atan(localPosition.z, localPosition.x);
        float u = (theta + PI) / (2.0 * PI);
        float v = 1.0 - clamp((localPosition.y + 1.2) / 2.4, 0.0, 1.0);
        surfaceData.surfaceColor = texture(iChannel0, vec2(u, v)).rgb;
    }


    // Se atingiu diretamente o sol, retorna apenas a cor intensa dele
    if (surfaceData.objectId == 30) {
        return sunColor * 1.0;
    }

    // Textura esférica para a Lua (ID 40)
    if (surfaceData.objectId == 40) {
        // Posição da lua como em getDist
        float sunAngle = iTime * 0.2;
        float orbitRadius = 6.0;
        vec3 moonPosition = vec3(-orbitRadius * cos(sunAngle), -3.0 * sin(sunAngle) + 1.0, 4.0);

        // Direção normalizada do ponto na lua em relação ao centro
        vec3 localPosition = normalize(surfacePosition - moonPosition);

        // Coordenadas UV com mapeamento esférico adequado
        float u = 0.5 + atan(localPosition.z, localPosition.x) / (2.0 * PI);
        float v = acos(localPosition.y) / PI;

        surfaceData.surfaceColor = texture(iChannel1, vec2(u, v)).rgb;

        // Brilho sutil adicional
        surfaceData.surfaceColor += vec3(0.1, 0.1, 0.2);
    }


    // Composição final da luz (ambiente + difusa + especular)
    vec3 finalColor =
        surfaceData.surfaceColor * surfaceData.ambientCoeff +
        (surfaceData.surfaceColor * sunColor * sunDiffuse * surfaceData.diffuseCoeff) * sunShadowFactor +
        (surfaceData.surfaceColor * redLightColor * redDiffuse * surfaceData.diffuseCoeff) * redShadowFactor +
        (surfaceData.surfaceColor * moonColor * moonDiffuse * surfaceData.diffuseCoeff * 0.4) + // lua mais fraca
        sunSpecular + redSpecular + moonSpecular * 0.4;


    return finalColor;
}

void main() {

    // Coordenadas normalizadas do pixel no plano de projeção da câmera
    vec2 screenUV = (gl_FragCoord.xy - 0.5 * iResolution) / iResolution.y;

    // Direção inicial do raio (antes da rotação pela câmera)
    vec3 rayDirection = normalize(vec3(screenUV, 1));

    // Cor final do pixel
    vec3 finalColor;

    // Coordenadas normalizadas do mouse (0~1)
    vec2 mouseUV = iMouse.xy / iResolution.xy;

    // Ângulos de rotação da câmera com base no mouse
    float azimuthAngle = mix(-0.5 * PI, 0.5 * PI, iMouse.z * mouseUV.x);    // horizontal
    float elevationAngle = mix(-0.25 * PI, 0.25 * PI, iMouse.z * mouseUV.y); // vertical

    // Ponto da caixa d'água (sempre visado)
    vec3 targetPosition = vec3(0.0, 1.2, 4.0);

    // Cálculo da altura do Sol e da Lua
    float sunAngle = iTime * 0.2;
    float sunY  =  3.0 * sin(sunAngle);
    float moonY = -3.0 * sin(sunAngle);

    // Determina o astro ativo acima do solo
    float activeY = (sunY > 0.0) ? sunY : (moonY > 0.0 ? moonY : 0.0);

    // Movimento inverso da câmera em relação ao astro ativo
    float cameraYOffset = -activeY * 0.5;
    float baseCameraY = 1.8;
    float cameraY = baseCameraY + cameraYOffset;

    // Raio da órbita da câmera ao redor do alvo
    float CameraOrbitRadius = 6.0;

    // Posição da câmera orbitando ao redor do ponto-alvo
    vec3 cameraPosition = targetPosition +
        vec3(cos(azimuthAngle), 0.0, sin(azimuthAngle)) * CameraOrbitRadius;

    // Aplica altura ajustada da câmera
    cameraPosition.y = cameraY;

    // Gera a matriz de orientação da câmera
    mat3 cameraMatrix = setCamera(cameraPosition, targetPosition);

    // Aplica a rotação da câmera à direção do raio
    rayDirection = cameraMatrix * rayDirection;

    // Executa ray marching a partir da câmera nessa direção
    Surface surfaceData = rayMarching(cameraPosition, rayDirection);

    // Ponto de interseção com a superfície (se encontrado)
    vec3 hitPosition = cameraPosition + surfaceData.distanceToSurface * rayDirection;

    // Calcula iluminação no ponto de impacto
    vec3 lighting = getLight(hitPosition, surfaceData, cameraPosition);

    // Define a cor do pixel com base na colisão ou fundo
    if (surfaceData.distanceToSurface < MAX_DIST) {
        finalColor = lighting; // Superfície visível
    } else {
        finalColor = surfaceData.surfaceColor; // Usa cor do céu com interpolação
    }

    // Correção gama para melhor aparência visual
    finalColor = pow(finalColor, vec3(0.4545));

    // Saída da cor final do fragmento
    C = vec4(finalColor, 1.0);
}
