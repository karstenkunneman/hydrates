Usage Examples
====================================================================================

The following code will return the equilibrium pressure (in Pa) of a pure methane hydrate system at 280K:

.. code-block:: python

    >> eqPressure = round(model.KlaudaSandler2003([1], [1], "T", 280, None).pressure)
    >> eqPressure
    >> 5145599