# Pruebas

```bash
conda activate cribado && pytest tests/ -q
```

Tarda unos segundos. No necesita Vina, ADCP, gnina ni conexión: cubre la lógica, no los motores.

## Qué cubre

- **`test_core.py`** — química de péptidos (ciclación cabeza-cola, protección de extremos, carga
  neta), traducción de rutas de Windows, reparto de recursos por memoria y flexibilidad,
  reconocimiento de péptidos por estructura y exportación en memoria.
- **`test_interfaz.py`** — dibuja cada etapa y cada modo de ligandos buscando excepciones, con la
  carpeta vacía y con datos, y comprueba que el estado sobrevive al cambiar de etapa.

Cada prueba documenta en su cadena de texto **el fallo real que la motivó**. Una prueba sin ese
contexto acaba borrándose cuando estorba.

## Qué NO cubre

Lo que depende de binarios externos y de una diana concreta: el acoplamiento en sí, PLIP, la
preparación con AGFR y la predicción ADMET. Para eso hace falta una corrida real; estas pruebas
solo garantizan que nada de lo que las rodea está roto.

## Al añadir una

Escribe primero la prueba que reproduce el fallo, compruébala en rojo y arregla después. Las tres
cuartas partes de los fallos de esta aplicación han sido silenciosos —un resultado plausible pero
equivocado, no una excepción—, así que comprueba el **valor**, no solo que no reviente.
