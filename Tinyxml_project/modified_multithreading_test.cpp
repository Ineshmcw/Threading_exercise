#include "./tinyxml2_modified/tinyxml2.h"
#include <iostream>
#include <vector>
// #include <chrono>
#include <unistd.h>

using namespace tinyxml2;

void manipulateXML(XMLDocument& doc, const std::string& inputFile) {
    XMLElement* root = doc.RootElement();
    if (root) {
        for (XMLElement* person = root->FirstChildElement("person"); person; person = person->NextSiblingElement("person")) {
            int age = 0;
            person->QueryIntAttribute("age", &age);
            person->SetAttribute("age", age + 1); // Increment age
        }
    }

    // Derive output file name
    std::string outputFile = "./output_files_modified/modified_" + inputFile.substr(inputFile.find_last_of("/\\") + 1);

    // Save the modified file
    if (doc.SaveFile(outputFile.c_str()) != XML_SUCCESS) {
        std::cerr << "Failed to save modified file: " << outputFile << std::endl;
    } else {
        // Update the print statement to show processed file paths
        std::cout << "Processed: " << inputFile << " -> " << outputFile << std::endl;
    }
}

int main() {
    std::vector<std::string> filePaths = {
        "./xml_files/example.xml",
        "./xml_files/file1.xml",
        "./xml_files/file2.xml"
    };

    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::cout << "Current working directory: " << cwd << std::endl;
    } else {
        std::cerr << "Error getting current working directory!" << std::endl;
    }

    // Use an index-based approach to pass file names contextually
    size_t index = 0;
    XMLDocument::LoadFiles(filePaths, [&filePaths, &index](XMLDocument& doc) mutable {
        std::string inputFile = filePaths[index++];
        manipulateXML(doc, inputFile);
    });

    return 0;
}
