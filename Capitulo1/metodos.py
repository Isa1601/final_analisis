import math
import numpy as np
import pandas as pd

def f_expr(expr, x):
    return eval(expr, {"x": x, "math": math})

def _compute_all_errors(prev, curr, dfx_str=None, cond_x=None):
    """
    Calcula y devuelve todos los errores:
      - 'abs'  : error absoluto = abs(curr - prev) (None si prev es None)
      - 'rel1' : abs((curr - prev) / curr)  (None si prev es None; inf si curr == 0)
      - 'rel2' : abs((curr - prev) / prev)  (None si prev es None; inf si prev == 0)
      - 'cond' : abs(f'(cond_x)) si dfx_str se proporciona; en caso contrario None

    prev, curr: aproximaciones anteriores y actuales (floats o None)
    dfx_str: expresión (string) de la derivada para evaluar el error de condición
    cond_x: punto donde evaluar la derivada (si None se usa curr)
    """
    abs_e = None
    rel1 = None
    rel2 = None
    cond = None

    if prev is not None and curr is not None:
        try:
            abs_e = abs(curr - prev)
        except Exception:
            abs_e = None

        try:
            if curr != 0:
                rel1 = abs((curr - prev) / curr)
            else:
                rel1 = float("inf")
        except Exception:
            rel1 = float("inf")

        try:
            if prev != 0:
                rel2 = abs((curr - prev) / prev)
            else:
                rel2 = float("inf")
        except Exception:
            rel2 = float("inf")

    # error de condición: requiere dfx_str
    if dfx_str is not None:
        x = cond_x if cond_x is not None else curr
        try:
            cond = abs(f_expr(dfx_str, x))
        except Exception:
            cond = None

    return {"abs": abs_e, "rel1": rel1, "rel2": rel2, "cond": cond}

def _min_non_null_error(err_dict):
    vals = [v for v in err_dict.values() if v is not None]
    return min(vals) if vals else None

def biseccion(fx_str, a, b, tol, niter, dfx_str=None):
    a, b, tol, niter = float(a), float(b), float(tol), int(niter)
    fa, fb = f_expr(fx_str, a), f_expr(fx_str, b)
    if fa * fb > 0:
        return {"error": "No hay cambio de signo en el intervalo [a, b]. Verifica que f(a) y f(b) tengan signos opuestos."}
    resultados = []
    prev_c = None
    for i in range(1, niter + 1):
        c = (a + b) / 2
        fc = f_expr(fx_str, c)

        # calcular todos los errores
        if prev_c is None:
            # primera iteración: no hay prev para relativo; usamos ancho del intervalo como aproximación del absoluto
            errores = {"abs": abs(b - a), "rel1": None, "rel2": None}
            if dfx_str is not None:
                try:
                    errores["cond"] = abs(f_expr(dfx_str, c))
                except Exception:
                    errores["cond"] = None
            else:
                errores["cond"] = None
        else:
            errores = _compute_all_errors(prev_c, c, dfx_str=dfx_str, cond_x=c)

        resultados.append({"iter": i, "a": a, "b": b, "c": c, "f(c)": fc, "errores": errores})

        # criterio de paro: si f(c) pequeño o cualquier error computado es menor que tol
        min_err = _min_non_null_error(errores)
        if abs(fc) < tol or (min_err is not None and min_err < tol):
            return {"resultados": resultados, "raiz": c}

        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
        prev_c = c

    return {"resultados": resultados, "error": "No converge dentro del número de iteraciones. Verifica el intervalo inicial y la tolerancia."}

def punto_fijo(gx_str, x0, tol, niter, dfx_str=None):
    x0, tol, niter = float(x0), float(tol), int(niter)
    resultados = []
    for i in range(niter):
        x1 = f_expr(gx_str, x0)
        errores = _compute_all_errors(x0, x1, dfx_str=dfx_str, cond_x=x1)
        resultados.append({"iter": i+1, "x": x0, "x1": x1, "errores": errores})
        min_err = _min_non_null_error(errores)
        if min_err is not None and min_err < tol:
            return {"resultados": resultados, "raiz": x1}
        x0 = x1
    return {"resultados": resultados, "error": "No converge dentro del número de iteraciones. Verifica que |g'(x)| < 1 cerca de la raíz y revisa el valor inicial."}

def regla_falsa(fx_str, a, b, tol, niter, dfx_str=None):
    a, b, tol, niter = float(a), float(b), float(tol), int(niter)
    fa, fb = f_expr(fx_str, a), f_expr(fx_str, b)
    if fa * fb > 0:
        return {"error": "No hay cambio de signo en el intervalo [a, b]. Verifica los extremos."}
    resultados = []
    prev_c = None
    for i in range(1, niter + 1):
        c = b - fb * (b - a) / (fb - fa)
        fc = f_expr(fx_str, c)

        if prev_c is None:
            # antes se usaba |f(c)| como "error" en la primera iteración; lo mantenemos para compatibilidad
            errores = {"abs": abs(fc), "rel1": None, "rel2": None}
            if dfx_str is not None:
                try:
                    errores["cond"] = abs(f_expr(dfx_str, c))
                except Exception:
                    errores["cond"] = None
            else:
                errores["cond"] = None
        else:
            errores = _compute_all_errors(prev_c, c, dfx_str=dfx_str, cond_x=c)

        resultados.append({"iter": i, "a": a, "b": b, "c": c, "f(c)": fc, "errores": errores})

        min_err = _min_non_null_error(errores)
        if abs(fc) < tol or (min_err is not None and min_err < tol):
            return {"resultados": resultados, "raiz": c}

        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
        prev_c = c
    return {"resultados": resultados, "error": "No converge dentro del número de iteraciones. Considera usar un mejor intervalo o menor tolerancia."}

def secante(fx_str, x0, x1, tol, niter, dfx_str=None):
    x0, x1, tol, niter = float(x0), float(x1), float(tol), int(niter)
    resultados = []
    for i in range(niter):
        f0 = f_expr(fx_str, x0)
        f1 = f_expr(fx_str, x1)
        if (f1 - f0) == 0:
            return {"error": "División por cero en la fórmula de la secante. Verifica que f(x0) ≠ f(x1)."}
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        errores = _compute_all_errors(x1, x2, dfx_str=dfx_str, cond_x=x2)
        resultados.append({"iter": i+1, "x0": x0, "x1": x1, "x2": x2, "errores": errores})
        min_err = _min_non_null_error(errores)
        if min_err is not None and min_err < tol:
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
        errores = _compute_all_errors(x0, x1, dfx_str=dfx_str, cond_x=x1)
        resultados.append({"iter": i+1, "x0": x0, "x1": x1, "f(x0)": f, "df(x0)": df, "errores": errores})
        min_err = _min_non_null_error(errores)
        if min_err is not None and min_err < tol:
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

        # método para raíces múltiples (ejemplo simple, m=1 equivalente a Newton modificado según contexto)
        x1 = x0 - (f / df)
        errores = _compute_all_errors(x0, x1, dfx_str=dfx_str, cond_x=x1)

        resultados.append({
            "iter": iteracion + 1,
            "x0": x0,
            "x1": x1,
            "f(x0)": f,
            "df(x0)": df,
            "ddf(x0)": ddf,
            "errores": errores
        })

        min_err = _min_non_null_error(errores)
        if abs(f) < tol or (min_err is not None and min_err < tol):
            return {"resultados": resultados, "raiz": x1}
        x0 = x1

    return {"resultados": resultados, "error": "No converge dentro del número de iteraciones. Verifica derivadas y valor inicial."}
