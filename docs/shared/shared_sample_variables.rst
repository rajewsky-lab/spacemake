One of the most important parts of spacemake are the so-called 'shared sample-variables'.
These are reusable, user-definable variables, which we can assign to several samples.
They can be shortly defined as follows:

``species``
   a collection of genome, annotation and rRNA\_genome. There is no default species, and each sample can have exactly one species.

``barcode-flavor``
   the variable which specifies the structure of Read1 and Read2, namely how the cell-barcode and UMI should be extracted. If no value provided for a sample, the default will be used.

``run-mode``
   each sample can have several ``run-mode``s, all of which are user definable. If no ``run-mode``s are specified, a sample will be processed using ``default`` ``run-mode`` settings.

``puck`` (spatial only)
   if a sample is spatial, it has to have a puck variable. If no puck is specified, a default puck will be used.  


To add, update, delete or list a shared sample-variable, you can use the following commands::

   spacemake config add-<shared-sample-variable>
   spacemake config update-<shared-sample-variable>
   spacemake config delete-<shared-sample-variable>
   spacemake config list-<shared-sample-variable>

where ``<shared-sample-variable>`` is one of ``species, barcode-flavor, run-mode or puck``
