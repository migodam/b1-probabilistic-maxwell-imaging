"""Standard-library runner; no installation into frozen A25 environment."""
import unittest,importlib.util,inspect,tempfile,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'src'))
suite=unittest.defaultTestLoader.discover(str(root/'tests'),pattern='test_analytic.py')
for name in ['test_contracts','test_vie_algebra']:
    spec=importlib.util.spec_from_file_location(name,root/'tests'/f'{name}.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    for function_name,function in vars(mod).items():
        if not function_name.startswith('test_') or not callable(function):continue
        def run(fn=function):
            with tempfile.TemporaryDirectory(prefix='research-w-test-') as tmp:
                if 'tmp_path' in inspect.signature(fn).parameters:fn(Path(tmp))
                else:fn()
        suite.addTest(unittest.FunctionTestCase(run,description=f'{name}.{function_name}'))
result=unittest.TextTestRunner(verbosity=2).run(suite)
sys.exit(0 if result.wasSuccessful() else 1)
