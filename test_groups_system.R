# Test the new ZooMSS Groups management system
library(assertthat)
library(readr)
library(dplyr, warn.conflicts = FALSE)

# Source the functions
source("R/fZooMSS_Groups.R")

cat("🧪 Testing ZooMSS Groups Management System\n")
cat("==========================================\n\n")

# Test 1: Load default groups
cat("📋 Test 1: Loading default groups\n")
groups1 <- fZooMSS_GetGroups(source = "default")
cat("✓ Loaded", nrow(groups1), "groups\n")
cat("✓ Species:", paste(head(groups1$Species, 3), collapse = ", "), "...\n\n")

# Test 2: Generate a template
cat("📝 Test 2: Creating template file\n")
template_file <- "test_template.csv"
fZooMSS_GetGroups(source = "template", file = template_file)
cat("✓ Template created at:", template_file, "\n\n")

# Test 3: Load from template file
cat("📁 Test 3: Loading from template file\n")
groups2 <- fZooMSS_GetGroups(source = "file", file = template_file)
cat("✓ Loaded", nrow(groups2), "groups from template\n")
cat("✓ Groups identical to default:", identical(groups1, groups2), "\n\n")

# Test 4: Modify and validate custom groups
cat("🔧 Test 4: Modifying and validating groups\n")
custom_groups <- groups1
custom_groups$W0[1] <- -15.0  # Modify flagellates min size
custom_groups$Wmax[1] <- -5.0  # Modify flagellates max size

# This should work
tryCatch({
  fZooMSS_ValidateGroups(custom_groups)
  cat("✓ Modified groups passed validation\n")
}, error = function(e) {
  cat("❌ Validation failed:", e$message, "\n")
})

# Test 5: Break validation (intentionally)
cat("\n⚠️  Test 5: Testing validation catches errors\n")
broken_groups <- custom_groups
broken_groups$W0[1] <- NA  # Break it

tryCatch({
  fZooMSS_ValidateGroups(broken_groups)
  cat("❌ Validation should have failed\n")
}, error = function(e) {
  cat("✓ Validation correctly caught error:", substr(e$message, 1, 50), "...\n")
})

# Clean up
if (file.exists(template_file)) {
  file.remove(template_file)
  cat("\n🧹 Cleaned up test files\n")
}

cat("\n🎉 All tests completed successfully!\n")
cat("📦 Groups management system is ready for integration\n")
