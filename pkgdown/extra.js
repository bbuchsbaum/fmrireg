(function () {
  var FAMILY_CLASSES = ["red", "lapis", "ochre", "teal", "green", "violet"];
  var PRESET_CLASSES = ["homage", "study", "structural", "adobe", "midnight"];
  var STYLE_CLASSES = ["minimal", "assertive"];

  function removeClasses(values, prefix) {
    values.forEach(function (value) {
      document.body.classList.remove(prefix + value);
    });
  }

  function applyDefaults() {
    if (!document.body) return;

    removeClasses(FAMILY_CLASSES, "palette-");
    removeClasses(PRESET_CLASSES, "preset-");
    removeClasses(STYLE_CLASSES, "style-");

    document.body.classList.add("palette-red", "preset-homage", "style-minimal");

    var theme = document.body.classList.contains("preset-midnight") ? "dark" : "light";
    document.documentElement.setAttribute("data-bs-theme", theme);
    document.body.setAttribute("data-bs-theme", theme);

    var nav = document.querySelector("nav.navbar");
    if (nav) nav.setAttribute("data-bs-theme", theme);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", applyDefaults);
  } else {
    applyDefaults();
  }
})();

