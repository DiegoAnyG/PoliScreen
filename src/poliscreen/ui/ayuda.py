"""Textos de ayuda de la interfaz.

Se mantienen fuera de streamlit_app.py por dos motivos: la interfaz queda legible, y las
explicaciones largas dejan de repetirse junto a cada control. En pantalla solo va una línea
breve; el desarrollo completo vive aquí y se consulta desde Ayuda.

Estructura: SECCIONES = {sección: [(titulo, cuerpo markdown), ...]}
"""

# Frases cortas que acompanan a un control concreto. La explicacion larga esta en SECCIONES.
BREVE = {
    "confianza": "Cuánto fiarse del resultado, no cuán bueno es. Detalle en Ayuda › Resultados.",
    "efectividad": "Porcentaje respecto al ligando de referencia del sitio. Detalle en Ayuda › Resultados.",
    "pki": "−log₁₀ Ki: numérico y ordenable. Detalle en Ayuda › Resultados.",
    "le": "Afinidad por átomo pesado. Detalle en Ayuda › Resultados.",
    "pose": "Todos los valores de la fila salen de esta misma pose.",
    "cavidades": "Extensión real frente a caja de búsqueda. Detalle en Ayuda › Ejecutar.",
    "peptidos": "Reglas y descriptores en Ayuda › Ligandos.",
    "gnina": "Segunda puntuación independiente. Detalle en Ayuda › Ejecutar.",
}

SECCIONES = {
    "Primeros pasos": [
        ("Qué hace PoliScreen",
         "Encadena el ciclo completo de cribado: preparar la diana, obtener compuestos, acoplarlos "
         "y evaluar la calidad de la unión.\n\n"
         "La diferencia con un panel de docking corriente está en cómo puntúa. En vez de premiar "
         "la afinidad y el número de contactos, mide **qué contactos** se hacen y **con qué residuos**, "
         "y añade una medida de **cuánto fiarse** de cada resultado."),
        ("El orden de trabajo",
         "1. **Receptores** — descarga o sube la estructura, elige qué conservar y extrae el ligando "
         "co-cristalizado como control.\n"
         "2. **Ligandos** — construye la serie por reacción, genera péptidos o sube estructuras listas.\n"
         "3. **Ejecutar** — define dónde buscar y lanza el acoplamiento.\n"
         "4. **Resultados** — ajusta la ponderación y examina el ranking.\n\n"
         "El panel derecho muestra siempre lo que corresponde a la etapa activa."),
        ("Guardar y recuperar el trabajo",
         "**Archivo › Guardar sesión** empaqueta el análisis en un archivo `.poliscreen`. Al restaurarlo "
         "vuelven las tablas, los receptores y los ligandos, y puedes **cambiar la ponderación sin "
         "repetir el docking**.\n\n"
         "La sesión ligera ocupa unos pocos megabytes; la completa añade poses y complejos."),
        ("La carpeta del proyecto",
         "Todo se escribe ahí: poses, complejos, XML de PLIP y tablas. Cambiar de carpeta cambia de "
         "análisis, y lo que ya haya preparado en ella (receptores, controles, ligandos) se detecta "
         "solo al entrar.\n\n"
         "PoliScreen se ejecuta dentro de Linux, así que la ruta es de Linux "
         "(`/home/usuario/...`). Si pegas una de Windows —`\\\\wsl.localhost\\...`, `C:\\...` o el "
         "trozo `home\\usuario\\...` que copia el Explorador— se traduce sola y se te avisa de cuál "
         "se ha usado."),
        ("Qué conviene descargar",
         "**Archivo › Descargar resultados** arma un solo ZIP con lo que marques, sin dejar copias en "
         "la carpeta del proyecto.\n\n"
         "Casi todo lo que aparece ahí **ya está en esa carpeta**: descargarlo solo tiene sentido para "
         "llevarte el análisis a otra máquina, adjuntarlo a un manuscrito o archivarlo. Por eso están "
         "marcados como **recomendados** los elementos que forman ese paquete mínimo:\n\n"
         "- La **sección de Métodos**, que es el único que no existe hasta que lo exportas.\n"
         "- El **ranking**, la **matriz de interacciones** y la **tabla de ligandos**, que son el "
         "resultado y su trazabilidad.\n"
         "- La **validación por redocking**, lo primero que revisa un evaluador.\n"
         "- **Receptores y ligandos de entrada**, para que otra persona pueda repetir la corrida.\n\n"
         "Lo demás —poses, complejos fusionados, resumen y energías de todas las poses— se regenera "
         "volviendo a ejecutar, y puede ocupar cientos de megabytes."),
    ],
    "Receptores": [
        ("Qué hace la preparación",
         "Elimina las aguas, añade hidrógenos y conserva la numeración original de los residuos.\n\n"
         "Conservar la numeración importa: si cambia, un residuo del sitio activo puede acabar "
         "identificado con el nombre de otro y todo el análisis farmacofórico queda mal."),
        ("Cofactores: cuándo conservarlos",
         "Un cofactor que forma parte del sitio (NADP en una reductasa, por ejemplo) debe conservarse: "
         "el ligando compite o coopera con él y quitarlo cambia la forma del bolsillo.\n\n"
         "Los iones y las moléculas de cristalización que no participan se pueden retirar."),
        ("El control co-cristalizado",
         "Es la pieza más importante del montaje. Define la referencia contra la que se mide todo, "
         "marca el sitio real de unión y permite validar el protocolo por redocking.\n\n"
         "Al extraerlo conviene indicar su **SMILES**: el formato PDB no guarda los órdenes de enlace y, "
         "sin esa plantilla, heterociclos como el N-óxido del benzofuroxano se leen mal."),
        ("Por qué se quitan las aguas",
         "Es la práctica habitual en acoplamiento rígido. AutoDock Vina no modela aguas explícitas ni "
         "el coste de desplazarlas, así que dejarlas produce artefactos: bloquean el sitio o generan "
         "contactos que no son reales.\n\n"
         "La excepción son las **aguas estructurales conservadas** que median la unión en algunos "
         "sistemas. Conservarlas es una decisión que hay que justificar caso por caso; si lo haces, "
         "PoliScreen detectará los puentes mediados por agua automáticamente."),
    ],
    "Ligandos": [
        ("Las tres vías",
         "**Construir por reacción** — parte de un núcleo y una biblioteca de reactivos, y filtra por "
         "viabilidad química real. Lo que se acopla es lo que se puede sintetizar.\n\n"
         "**Generar péptidos** — enumera secuencias bajo reglas de composición y propiedades.\n\n"
         "**Subir ligandos listos** — estructuras ya preparadas. Se lee su estructura química para "
         "poder calcular ADMET y descriptores."),
        ("Péptidos: las reglas",
         "El alfabeto se restringe por clases (hidrofóbicos, catiónicos, aromáticos…) y se pueden fijar "
         "un prefijo o un sufijo, prohibir repeticiones o limitar residuos consecutivos.\n\n"
         "Los filtros fisicoquímicos son los que más discriminan en péptidos antimicrobianos: la "
         "**carga neta positiva** (la membrana bacteriana es aniónica) y una hidrofobicidad moderada."),
        ("Péptidos: los descriptores",
         "**Carga neta** a pH 7.4 — el rasgo más asociado a la actividad antimicrobiana.\n\n"
         "**GRAVY** — hidropatía media; positivo indica carácter hidrofóbico global.\n\n"
         "**Momento hidrofóbico** — mide la anfipaticidad: si al plegarse en hélice los residuos "
         "hidrofóbicos quedan en una cara y los polares en la otra. Es lo que permite insertarse en "
         "la membrana.\n\n"
         "**Índice de Boman** — tendencia a unirse a otras proteínas; por encima de 2.5 kcal/mol se "
         "considera promiscuo."),
        ("Péptidos: los extremos",
         "**Amidar el extremo C** elimina la carga negativa terminal y suma +1 a la carga neta, lo que "
         "suele aumentar la actividad antimicrobiana.\n\n"
         "**Acetilar el extremo N** protege frente a aminopeptidasas.\n\n"
         "**Ciclar cabeza-cola** rigidiza el péptido: reduce mucho los grados de libertad, resiste "
         "proteasas y además hace más fiable el acoplamiento."),
        ("Límite del acoplamiento de péptidos",
         "Medido sobre saFtsZ con AutoDock Vina (caja de 23 Å, exhaustividad 8, un hilo):\n\n"
         "| Residuos | Enlaces rotables | Tiempo |\n|---|---|---|\n"
         "| 3 | 15 | ~98 s |\n| 5 | 23 | más de 2 min |\n| 10 | 43 | no termina |\n\n"
         "Es una limitación de Vina, no de PoliScreen: trata el ligando como un árbol de torsiones "
         "independientes y con muchas el muestreo deja de cubrir el espacio. **Para péptidos, el "
         "acoplamiento sirve para ordenar candidatos, no para proponer un modo de unión.**"),
    ],
    "Ejecutar": [
        ("La caja de búsqueda",
         "Lo más fiable es centrarla en el ligando co-cristalizado: marca el sitio real. El centro "
         "geométrico de la proteína o un cofactor apuntan a otro lugar.\n\n"
         "Los ejes X, Y y Z se dibujan en el visor para saber en qué dirección mueve cada control."),
        ("Cavidades detectadas",
         "`Cavidad` es la extensión real del bolsillo; `Caja` es la región de búsqueda que se le asigna, "
         "con un **mínimo de 14 Å** porque por debajo no cabría un ligando. Cuando se aplica ese mínimo "
         "se marca con `*`, y por eso cavidades de volumen distinto pueden compartir el mismo tamaño de caja.\n\n"
         "La **drogabilidad** estima si el bolsillo tiene forma y química adecuadas para unir una molécula "
         "pequeña. `Flexibility` no se muestra porque se deriva de los factores B, que la preparación deja a cero."),
        ("Docking híbrido",
         "Acopla los mismos compuestos en **varios bolsillos** del mismo receptor, cada uno con su ranking.\n\n"
         "Sirve para responder si un compuesto prefiere el sitio catalítico o se cuela en uno alostérico: "
         "es información de **selectividad de sitio** que un cribado de un solo sitio no da.\n\n"
         "Cada sitio usa su propia referencia: el control si está ahí, un cofactor si cae dentro de la caja, "
         "o los residuos catalíticos que designes."),
        ("Parámetros de acoplamiento",
         "**Exhaustividad** — cuánto explora la búsqueda. Más alto es más fino y más lento.\n\n"
         "**Poses por ligando** — cuántos modos de unión se conservan. Por debajo de 3 la métrica de "
         "confianza pierde resolución.\n\n"
         "**Rango de energía** — ventana respecto a la mejor pose para reportar modos alternativos. No "
         "cambia la mejor pose ni la profundidad de búsqueda.\n\n"
         "**Semilla y un hilo** — garantizan que dos corridas iguales den el mismo resultado. Con más de "
         "un hilo por acoplamiento, Vina deja de ser determinista."),
        ("Segunda puntuación con gnina",
         "Vuelve a evaluar con una **red neuronal** las poses que Vina ya generó, sin repetir la búsqueda.\n\n"
         "Su valor no es la velocidad sino la independencia: la puntuación de Vina es empírica (términos "
         "físicos ajustados a datos experimentales) y la de gnina se aprende de complejos cristalográficos. "
         "Que dos métodos con supuestos distintos coincidan es una evidencia que ninguno da por separado, "
         "y esa concordancia entra en la métrica de confianza.\n\n"
         "**Limitación**: la red solo puede juzgar las poses que Vina encontró. Si el muestreo no dio con "
         "la pose correcta, el re-puntuado no la recupera."),
    ],
    "Resultados": [
        ("Efectividad",
         "Porcentaje respecto al ligando de referencia de ese sitio, que queda en 100 %.\n\n"
         "La calidad de interacción **no mide parecido al control**: suma cada contacto ponderado por su "
         "tipo de enlace (salino > puente de hidrógeno > π > hidrofóbico) y por el rol del residuo "
         "(catalítico, secundario, de cavidad o externo). Así se supera al control haciendo **más y "
         "mejores contactos productivos**, no copiándolo."),
        ("Confianza: qué es y por qué es distinta",
         "Es **ortogonal a la efectividad**: no mide cuán bueno es el compuesto, sino cuánto fiarse del número.\n\n"
         "Se calcula como media **geométrica** —basta que una evidencia falle para que baje— de:\n\n"
         "- **conv**: convergencia del modo de unión, el solape de interacciones entre las mejores poses.\n"
         "- **conc**: concordancia entre afinidad e interacción.\n"
         "- **consenso**: acuerdo entre Vina y la red neuronal, si se activó el re-puntuado.\n\n"
         "Se reduce si la diana no valida su redocking. **Un compuesto con efectividad alta y confianza "
         "baja es una alarma**: buen score, pero evidencias que no concuerdan.\n\n"
         "`geom` (dispersión geométrica de las poses) se muestra como diagnóstico pero no entra en el "
         "cálculo: con Vina resulta casi constante y no discrimina."),
        ("Las métricas de afinidad",
         "**best_dock** — energía de Vina en kcal/mol; más negativo es mejor.\n\n"
         "**pKi** — −log₁₀ de la Ki estimada. Es numérico y ordenable, a diferencia de la columna Ki, que "
         "es texto con unidades mezcladas. La Ki **no entra en el puntaje**: se deriva del score de "
         "docking y puntuarla sería contar lo mismo dos veces.\n\n"
         "**LE** (−ΔG/átomos pesados) y **LLE** (pKi−LogP) premian la unión **por átomo**. Sirven de guarda "
         "contra el sesgo de tamaño: sin ellas, la molécula más grande casi siempre gana."),
        ("Residuos catalíticos y secundarios",
         "**Catalíticos** — obligatorios: no tocarlos penaliza mediante la puerta catalítica. Son una "
         "propiedad de la enzima, determinada experimentalmente, no del ligando.\n\n"
         "**Secundarios** — anclas conocidas del bolsillo que suman más que un contacto cualquiera, pero "
         "no se exigen.\n\n"
         "Si no conoces el sitio catalítico, PoliScreen sugiere los residuos con los que el ligando "
         "cristalográfico hace interacciones direccionales. Es un punto de partida, no una respuesta: "
         "**esta lista influye en el ranking más que cualquier peso**."),
        ("Ponderación",
         "Los pesos de eje se **auto-normalizan**: no hace falta que sumen 1, y ponerlos todos a 1.0 "
         "equivale a un promedio simple. Si los pones todos a 0 no hay puntaje y la aplicación avisa.\n\n"
         "Un eje con peso pero **sin datos** (toxicidad sin haber predicho ADMET) se ignora, y también se avisa: "
         "así la sección de Métodos no declara algo que no intervino."),
        ("La columna pose",
         "Indica de qué modelo salen **todos** los valores de esa fila. Docking, interacciones y Ki vienen "
         "de la misma pose, nunca de modelos distintos."),
        ("Percentil frente a porcentaje",
         "El **porcentaje** se mide contra el control de esa diana, así que depende de lo fuerte que sea ese "
         "control: no es comparable entre dianas.\n\n"
         "El **percentil** sitúa al compuesto dentro de su propia biblioteca y sí permite comparar."),
    ],
    "Reproducibilidad": [
        ("Qué se fija",
         "Semilla constante y un hilo por acoplamiento; versiones fijadas del entorno y del binario de "
         "AutoDock Vina, verificado por su suma SHA256.\n\n"
         "Dos corridas con la misma configuración dan el mismo resultado."),
        ("Exportar los métodos",
         "**Archivo › Exportar Métodos** genera un documento con los parámetros, la caja, los pesos, la "
         "referencia empleada y las versiones exactas de cada herramienta, listo para la sección de "
         "Métodos de un artículo."),
        ("Limitaciones declaradas",
         "Acoplamiento **rígido**: el receptor no cambia de conformación.\n\n"
         "**No covalente**: no se modelan enlaces covalentes con la diana.\n\n"
         "**Ki estimada** a partir de la afinidad de Vina, que es una aproximación de ΔG: es informativa, "
         "no medida.\n\n"
         "**LD50 y toxicidad** proceden de un solo modelo predictivo y tienden a ser optimistas. Conviene "
         "contrastarlas con un segundo predictor antes de afirmar baja toxicidad."),
    ],
}
