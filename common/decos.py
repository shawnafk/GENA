"""Small reusable decorators used across the GENA codebase."""

from __future__ import annotations

import functools
import time


def timer(func):
    """Print function execution time in seconds."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        run_time = time.perf_counter() - start_time
        print(f"Function '{func.__name__}' executed,耗时: {run_time:.4f} s")
        return result
    return wrapper


def logger(func):
    """Log function entry and exit."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        print(f"Preparing to execute: {func.__name__}")
        result = func(*args, **kwargs)
        print(f"Function {func.__name__} executed.")
        return result
    return wrapper


def repeat(num_times):
    """Return a decorator that repeats a function call ``num_times`` times."""
    def decorator_repeat(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            for _ in range(num_times):
                result = func(*args, **kwargs)
            return result
        return wrapper
    return decorator_repeat
