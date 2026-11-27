// On the articles index, wrap each article title and description in a card
// so the articles can be shown as a grid.
document.addEventListener("DOMContentLoaded", function () {
  if (!/\/articles\/(index\.html)?$/.test(window.location.pathname)) return;
  document.body.classList.add("bison-articles-index");
  document.querySelectorAll("main dl").forEach(function (dl) {
    var terms = Array.prototype.slice.call(dl.querySelectorAll(":scope > dt"));
    terms.forEach(function (dt) {
      var card = document.createElement("div");
      card.className = "bison-article-card";
      var dd = dt.nextElementSibling;
      dl.insertBefore(card, dt);
      card.appendChild(dt);
      if (dd && dd.tagName === "DD") card.appendChild(dd);
    });
    dl.classList.add("bison-article-grid");
  });
});
