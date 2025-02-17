plugins {
    kotlin("jvm") version "2.1.0"
    application
}

group = "org.example"
version = "1.0-SNAPSHOT"

repositories {
    mavenCentral()
}

dependencies {
    val javaFxVersion = "23.0.2"
    val platform = "mac-aarch64"  // Use ARM64 version

    implementation("org.openjfx:javafx-base:$javaFxVersion:$platform")
    implementation("org.openjfx:javafx-controls:$javaFxVersion:$platform")
    implementation("org.openjfx:javafx-graphics:$javaFxVersion:$platform")
    // Uncomment if you need FXML support:
    // implementation("org.openjfx:javafx-fxml:$javaFxVersion:$platform")
}

tasks.test {
    useJUnitPlatform()
}

application {
    mainClass.set("GraphTetrahedronAppKt")
    applicationDefaultJvmArgs = listOf(
        "--module-path", "/Users/ayasamadzelkava/Applications/javafx-sdk-23.0.2/lib",
        "--add-modules", "javafx.controls,javafx.fxml,javafx.graphics"
        // Do not force software rendering here.
    )
}

kotlin {
    jvmToolchain(21)
}
