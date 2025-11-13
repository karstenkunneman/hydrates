Usage Examples
====================================================================================

The following code will return the equilibrium pressure (in Pa) of a pure methane hydrate system at 280K:

.. code-block:: python

    >> eqPressure = round(p2f_HydrateCalcLib.model.KlaudaSandler2003([1], [1], "T", 280, None).pressure)
    >> eqPressure
    >> 5145599

The following code will return the equilibrium temperature (in K) of a methane-ethane-propane hydrate system at 10.7 bar:

.. code-block:: python

    >> eqTemperature = round(p2f_HydrateCalcLib.model.KlaudaSandler2003([1, 2, 3], [0.9196, 0.0513, 0.0291], "P", None, 10.7E5).temperature, 1)
    >> eqTemperature
    >> 275.8

