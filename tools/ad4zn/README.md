# AutoDock4Zn resources

This folder contains the supplied `zinc_pseudo.py`, `prepare_gpf4zn.py`, and
`AD4Zn.dat` files used by the optional AutoDock4Zn protocol. The application
discovers these files automatically from this checkout. Existing source and
copyright notices are retained in the supplied files.

ADFR Suite's compatible Python runtime (with MolKit and AutoDockTools) and
AutoGrid >= 4.2.7 must still be installed separately. Run the helper scripts
through that runtime; they are not application-environment dependencies.

See [setup and validation](../../docs/AD4Zn.md) for configuration, upstream
sources, and the required zinc docking validation. Standard Vina does not need
these resources.
