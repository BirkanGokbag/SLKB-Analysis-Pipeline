## Main API

### Sample Data

![Sample Data](SampleData.png)

Note: Label can be a number or be empty for non-classified celltypes. 

### preprocessing

Removes a set of columns if specified, returns input CYTOF data with labels, the labels as a list, and the unlabeled input CyTOF data.

```
X_data_labeled, y_data, data_unlabeled = DGCyTOF.preprocessing(dataset, columns_to_remove = [])
```

**Params**:

* dataset: CyTOF Data Matrix
* columns_to_remove: List of columns to remove from the dataset
    * default: Empty list

**Returns**:

* X_data_labeled: input CYTOF data with labels
* y_data: the labels as a list
* data_unlabeled: non-classified CyTOF data

<hr>

### train_model

Trains the entered deep learning model using Pytorch. Utilized criterion is CrossEntropyLoss and optimizer is Adam optimizer with learning rate 0.001. The input model is set to evaluation mode after training.

```
DGCyTOF.train_model(model_fc, X_train, max_epochs = 20, params_train = {'batch_size':128, 'shuffle': True, 'num_workers': 6})
```

**Params**:

* model_fc: PyTorch model, must have a forward function and utilize argmax as classification in its design
* max_epochs: Number of epochs the model will be trained for 
    * default: 20
* params_train: dictionary containing information for trainloader, requires at least a batch_size, shuffle, and num_workers keys.
    * batch_size: Number of data points in a single batch, default 128
    * shuffle: Shuffle the batches prior to training, default True
    * num_workers: Number of processes that will be used to load data, default 6

**Returns**:

* None: Nothing is returned, the model is set to eval mode after training. 

<hr>

### validate_model

Runs validation on the validation dataset, print out the performance of the trained model for all cell types and returns them as a zip.

```
validation_results = DGCyTOF.validate_model(model_fc, val_tensor, classes, params_val = {'batch_size':10000, 'shuffle': False, 'num_workers': 6})
```

**Params**:

* model_fc: Trained PyTorch model, must have a forward function and utilize argmax as classification in its design
* val_tensor: Validation dataset as a Torch tensor.
* classes: List of types of cells
* params_val: dictionary containing information for dataloader, requires at least a batch_size, shuffle, and num_workers keys.
    * batch_size: Number of data points in a single batch, default 128
    * shuffle: Shuffle the batches prior to training, default True
    * num_workers: Number of processes that will be used to load data, default 6

**Returns**:

* Zip of listed results. Each respective row contains pred,label,out in validation_results
    * pred: Predicted label of a data point
    * label: Actual label of a data point
    * out: Output value of the data running forward through model_fc

<hr>

### calibrate_data

Calibrates the test set of the data based on the training model and validation results in performing over the test set. Test data either are classified more accurately or labeled as a new subtype based on the minimum threshold in classification. Minumum threshold is computed by obtaining the lowest correlation probability from validation results.
Returns calibrated incorrect data.

```
updated_incorrect_data = DGCyTOF.calibrate_data(model_fc, X_test, classes, validation_results, unlabeled_data)
```

**Params**:

* model_fc: Trained PyTorch model, must have a forward function and utilize argmax as classification in its design
* X_test: CyTOF data for test set.
* classes: List of types of cells
* validation_results: zip created by validate_model. Contains the following keys:
    * pred: Predicted label of a data point
    * label: Actual label of a data point
    * out: Output value of the data running forward through model_fc

**Returns**:

* updated_incorrect_data: Calibrated incorrect data.

<hr>

### dimensionality_reduction_and_clustering

Reduce the dimensions of the input data using HDBSCAN + UMAP. While running, also displays the 2D clustered data (returned on clusterPlot).

```
data_umap, y_HDBSCAN_umap, new_subtypes, clusterPlot = DGCyTOF.dimensionality_reduction_and_clustering(input_data, n_neighbors = 5, min_dist = 0.01)
```

**Params**:

* input_data: CyTOF matrix data
* n_neighbors: Number of neighbors set to cluster data points
    * Default: 5
* min_dist: Minimum distance between embeded points, utilized in the UMAP function.
    * Default: 0.01

**Returns**:

* data_umap: Data with UMAP transform to 2 dimensions.
* y_HDBSCAN_umap: HDBSCAN + UMAP of the input data
* new_subtypes: Returns found new subtypes
* clusterPlot: 2D plot of data with new subtypes

<hr>

**Additional Functions**

The following two visualization functions are included in the package to allow the user to visualize their results quickly on the spot.

### Dim_Red_Plot

Creates a 2D plot of the data, displays it and returns it.

```
figure = DGCyTOF.Dim_Red_Plot(name, data, labels, no_classes, class_names)
```

**Params**:

* name: Name of the plot 
* data: 2 dimensional embedding of data, for instance UMAP transformed data
* no_classes: List of classes
* class_names: List of class names

**Returns**:

* figure: Plt plot

<hr>

### Dim_Red_Plot_3d

Creates a 3D plot of the data, displays it and returns it.

```
figure = DGCyTOF.Dim_Red_Plot_3d(data, labels, all_celltypes)
```

**Params**:

* data: Low dimension embedding of data, for instance UMAP transformed data
* labels: List of labels corresponding to each cell
* all_celltypes: List of celltype names, where each celltype is in its respective label position

**Returns**:

* figure: Plt plot

<hr>

© Copyright 2020, The Ohio State University, Biomedical Informatics