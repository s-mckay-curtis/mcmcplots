.html.end <- function(file) {
    script <- '
<script>
function filterPlots() {
  var input = document.getElementById("param_search");
  if (!input) return;
  var filter = input.value.trim().toLowerCase();
  var cards = document.querySelectorAll(".plot-card");
  var groups = document.querySelectorAll(".group-section");
  var tocItems = document.querySelectorAll(".toc_item");

  cards.forEach(function(card) {
    var param = (card.getAttribute("data-param") || "").toLowerCase();
    if (param.indexOf(filter) > -1) {
      card.style.display = "";
    } else {
      card.style.display = "none";
    }
  });

  groups.forEach(function(grp) {
    var visibleCards = grp.querySelectorAll(".plot-card:not([style*=\'display: none\'])");
    grp.style.display = visibleCards.length > 0 ? "" : "none";
  });

  tocItems.forEach(function(item) {
    var group = item.getAttribute("data-group");
    var grpElem = document.getElementById("group-" + group);
    if (grpElem && grpElem.style.display === "none") {
      item.style.display = "none";
    } else {
      item.style.display = "";
    }
  });
}
</script>
</body>
</html>
'
    cat(script, "\n", file = file, append = TRUE)
}
