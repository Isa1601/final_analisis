import math
import numpy as np
import pandas as pd

def f_expr(expr, x):
    return eval(expr, {"x": x, "math": math})

def biseccion(fx_str, a, b, tol, niter):
    a, b, tol, niter = float(a), float(b), float(tol), int(niter)
    fa, fb = f_expr(fx_str, a), f_expr(fx_str, b)
    if fa * fb > 0:
        return {"error": "No hay cambio de signo en el intervalo [a, b]. Verifica que f(a) y f(b) tengan signos opuestos."}
    resultados = []
    for i in range(1, niter + 1):
        c = (a + b) / 2
        fc = f_expr(fx_str, c)
        resultados.append({"iter": i, "a": a, "b": b, "c": c, "f(c)": fc, "error": abs(b - a)})
        if abs(fc) < tol or abs(b - a) < tol:
            return {"resultados": resultados, "raiz": c}
        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    return {"resultados": resultados, "error": "No converge dentro del número de iteraciones. Verifica el intervalo inicial y la tolerancia."}

def punto_fijo(gx_str, x0, tol, niter):
    x0, tol, niter = float(x0), float(tol), int(niter)
    resultados = []
    for i in range(niter):
        x1 = f_expr(gx_str, x0)
        err = abs(x1 - x0)
        resultados.append({"iter": i+1, "x": x0, "x1": x1, "error": err})
        if err < tol:
            return {"resultados": resultados, "raiz": x1}
        x0 = x1
    return {"resultados": resultados, "error": "No converge dentro del número de iteraciones. Verifica que |g'(x)| < 1 cerca de la raíz y revisa el valor inicial."}

def regla_falsa(fx_str, a, b, tol, niter):
    a, b, tol, niter = float(a), float(b), float(tol), int(niter)
    fa, fb = f_expr(fx_str, a), f_expr(fx_str, b)
    if fa * fb > 0:
        return {"error": "No hay cambio de signo en el intervalo [a, b]. Verifica los extremos."}
    resultados = []
    for i in range(1, niter + 1):
        c = b - fb * (b - a) / (fb - fa)
        fc = f_expr(fx_str, c)
        resultados.append({"iter": i, "a": a, "b": b, "c": c, "f(c)": fc, "error": abs(fc)})
        if abs(fc) < tol:
            return {"resultados": resultados, "raiz": c}
        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    return {"resultados": resultados, "error": "No converge dentro del número de iteraciones. Considera usar un mejor intervalo o menor tolerancia."}

def secante(fx_str, x0, x1, tol, niter):
    x0, x1, tol, niter = float(x0), float(x1), float(tol), int(niter)
    resultados = []
    for i in range(niter):
        f0 = f_expr(fx_str, x0)
        f1 = f_expr(fx_str, x1)
        if (f1 - f0) == 0:
            return {"error": "División por cero en la fórmula de la secante. Verifica que f(x0) ≠ f(x1)."}
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        err = abs(x2 - x1)
        resultados.append({"iter": i+1, "x0": x0, "x1": x1, "x2": x2, "error": err})
        if err < tol:
            return {"resultados": resultados, "raiz": x2}
        x0, x1 = x1, x2
    return {"resultados": resultados, "error": "No converge dentro del número de iteraciones. Prueba con x0 y x1 más cercanos a la raíz."}

def newton(fx_str, dfx_str, x0, tol, niter):
    x0, tol, niter = float(x0), float(tol), int(niter)
    resultados = []
    for i in range(niter):
        f = f_expr(fx_str, x0)
        df = f_expr(dfx_str, x0)
        if df == 0:
            return {"error": "Derivada cero, no se puede continuar. Verifica f'(x) cerca de la raíz."}
        x1 = x0 - f / df
        err = abs(x1 - x0)
        resultados.append({"iter": i+1, "x0": x0, "x1": x1, "f(x0)": f, "df(x0)": df, "error": err})
        if err < tol:
            return {"resultados": resultados, "raiz": x1}
        x0 = x1
    return {"resultados": resultados, "error": "No converge dentro del número de iteraciones. Verifica la derivada y el punto inicial."}

def raices_multiples(fx_str, dfx_str, ddfx_str, x0, tol, niter):
    x0, tol, niter = float(x0), float(tol), int(niter)
    resultados = []

    for iteracion in range(niter):
        try:
            f   = f_expr(fx_str,  x0)
            df  = f_expr(dfx_str, x0)
            ddf = f_expr(ddfx_str, x0)
        except Exception as e:
            return {"error": f"Error al evaluar funciones: {e}"}

        if df == 0:
            return {"error": f"Derivada cero en iteración {iteracion+1}. Verifica f'(x)."}

        x1 = x0 - (f / df)  # m = 1 en este caso
        err = abs(x1 - x0) #Error absoluto

        err1 = abs((x1-x0)/x1) #Error relativo 1

        err2 = abs((x1-x0)/x0) #Error relativo 2
        
        err3 = abs(dfx_str)
        
        resultados.append({
            "iter": iteracion + 1,
            "x0": x0,
            "x1": x1,
            "error abs": err
            "error rel1": err1
            "error rel2": err2
            "error con": err3
        })

        if abs(f) < tol or err < tol:
            return {"resultados": resultados, "raiz": x1}
        x0 = x1

    return {"resultados": resultados, "error abs": "No converge dentro del número de iteraciones. Verifica derivadas y valor inicial."}
