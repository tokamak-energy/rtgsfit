from rtgsfit_verify_analytic.analytic_soln import add

def test_add():
    assert add(2, 3) == 5
    assert add(-1, 1) == 0
