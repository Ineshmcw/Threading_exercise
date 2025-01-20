#include "./tinyxml2_original/tinyxml2.h"
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <filesystem> // Requires C++17
// #include <chrono> // For measuring execution time

using namespace tinyxml2;
namespace fs = std::filesystem;

void processXML(const std::string& inputPath, const std::string& outputPath) {
    XMLDocument doc;

    // Load the XML file
    if (doc.LoadFile(inputPath.c_str()) != XML_SUCCESS) {
        std::cerr << "Failed to load file: " << inputPath << " | " << doc.ErrorStr() << std::endl;
        return;
    }

    // Access the root element
    XMLElement* root = doc.RootElement();
    if (!root) {
        std::cerr << "No root element found in: " << inputPath << std::endl;
        return;
    }

    // Iterate over "person" elements and modify their "age"
    for (XMLElement* person = root->FirstChildElement("person"); person; person = person->NextSiblingElement("person")) {
        int age = 0;
        if (person->QueryIntAttribute("age", &age) == XML_SUCCESS) {
            person->SetAttribute("age", age + 1); // Increment age by 1
        }
    }

    // Save the modified XML to the output path
    if (doc.SaveFile(outputPath.c_str()) != XML_SUCCESS) {
        std::cerr << "Failed to save file: " << outputPath << " | " << doc.ErrorStr() << std::endl;
        return;
    }

    // Print output in the requested format
    std::cout << "Processed: " << inputPath << " -> " << outputPath << std::endl;
}

int main() {
    // auto startTime = std::chrono::high_resolution_clock::now(); // Start the timer

    std::string inputDir = "./xml_files/";
    std::string outputDir = "./output_files_original/";

    // Create output directory if it doesn't exist
    if (!fs::exists(outputDir)) {
        fs::create_directory(outputDir);
    }

    std::vector<std::thread> threads;

    // Iterate over XML files in the input directory
    for (const auto& entry : fs::directory_iterator(inputDir)) {
        if (entry.path().extension() == ".xml") {
            std::string inputPath = entry.path().string();
            // Modify the output filename by adding the "modified_" prefix
            std::string outputPath = outputDir + "modified_" + entry.path().filename().string();

            // Launch a thread to process each XML file
            threads.emplace_back(processXML, inputPath, outputPath);
        }
    }

    // Wait for all threads to complete
    for (auto& thread : threads) {
        thread.join();
    }

    // auto endTime = std::chrono::high_resolution_clock::now(); // End the timer
    // std::chrono::duration<double> elapsed_time = endTime - startTime;

    // std::cout << "All files processed in " << elapsed_time.count() << "seconds!" << std::endl;

    return 0;
}
