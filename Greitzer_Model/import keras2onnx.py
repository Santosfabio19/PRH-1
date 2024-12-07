import keras2onnx
keras2onnx.convert_keras(model, name=None, doc_string='', target_opset=None, channel_first_inputs=None)
    # type: (keras.Model, str, str, int, []) -> onnx.ModelProto
