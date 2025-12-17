from DataProcessor import DataProcessor
from Model import Model


if __name__ == "__main__":
    processor = DataProcessor("standard_lek")
    all_label_data = processor.extractFromFiles("out/")
    all_input_data = processor.generateInputData()
    data_scaled, input_scaled = processor.normalizeData(all_label_data, all_input_data)
    processor.train_test_val_split(data_scaled, input_scaled)
    # processor.printDataInfo()

    model = Model()
    model.load_data("standard_lek")
    model.initialize_model()
    model.train_model()
    model.save_model("standard_lek_model.pt")
    model.plot_loss()
