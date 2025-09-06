let slideIndex = 0;
const slides = document.querySelectorAll(".slide");
const dots = document.querySelectorAll(".dot");
const prevBtn = document.querySelector(".prev");
const nextBtn = document.querySelector(".next");

function showSlide(n) {
  slides.forEach((slide, i) => {
    slide.classList.remove("active");
    dots[i].classList.remove("active");
    if (i === n) {
      slide.classList.add("active");
      dots[i].classList.add("active");
    }
  });
  slideIndex = n;
}

prevBtn.addEventListener("click", () => {
  slideIndex = (slideIndex - 1 + slides.length) % slides.length;
  showSlide(slideIndex);
});

nextBtn.addEventListener("click", () => {
  slideIndex = (slideIndex + 1) % slides.length;
  showSlide(slideIndex);
});

dots.forEach((dot, i) => {
  dot.addEventListener("click", () => {
    showSlide(i);
  });
});

// Initialize
showSlide(slideIndex);
